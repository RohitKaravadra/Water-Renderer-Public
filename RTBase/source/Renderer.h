#pragma once

#include <thread>
#include <functional>

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Scene.h"
#include "GamesEngineeringBase.h"

#include "Denoiser.h"
#include "Settings.h"
#include "MediumData.h"

static const Color specrals[3]{
	Color(1.0f, 0.0f, 0.0f),	// Red
	Color(0.0f, 1.0f, 0.0f),	// Green
	Color(0.0f, 0.0f, 1.0f)		// Blue
};

class Renderer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom** samplers;					// samplers for multithreading
	std::thread** threads;					// threads for multithreading
	int numProcs;							// number of processors
	unsigned int numThreads;				// number of threads

	unsigned int totalTiles;				// number of tiles
	unsigned int totalXTiles;				// number of tiles in x direction	
	const unsigned int tileSize = 16;		// size of each tile

	std::atomic<unsigned int> tileCounter;	// number of tiles processed
	std::vector<unsigned int> tileSamples;	// number of samples per tile

	SETTINGS settings;						// settings for the ray tracer

	MediumData vlmData;				// water properties for volume rendering

	// structure to hold path data for path tracing
	struct PathData
	{
		Sampler* sampler;			// pointer to sampler

		Ray& r;						// Path ray
		Color pathThroughput;		// pathThroughput

		int depth;					// path depth
		float prevPdf;			// previous bsdf pdf;
		bool prevSpecular;			// previous bsdf specular

		// scattering event data
		float eta;											// Relative refractive index
		Color sigma_t;										// Extinction coefficient
		Color sigma_s;										// Scattering coefficient
		std::unique_ptr<PhaseFunction> scatterPhaseFunc;	// Owned phase function

		PathData(Ray& r, Sampler* sampler) :r(r), sampler(sampler)
		{
			pathThroughput = Color(1.0f, 1.0f, 1.0f);

			depth = 0;
			prevPdf = 0;
			prevSpecular = true;
		}

		~PathData() {
			sampler = nullptr;
		}
	};

	~Renderer()
	{
		std::cout << "Cleaning Ray Tracer..." << std::endl;

		// clean threads
		if (threads != nullptr)
			delete[] threads;

		// clean samplers
		if (samplers != nullptr)
		{
			for (unsigned int i = 0; i < numThreads; i++)
				delete samplers[i];
			delete[] samplers;
		}

		// clean film
		if (film != nullptr)
			delete film;

	}

	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas, SETTINGS _settings)
	{
		scene = _scene;
		canvas = _canvas;
		settings = _settings;

		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, settings.filter);
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;

		if (!(settings.adaptiveSampling && settings.useMultithreading))
			settings.initSPP = settings.totalSPP;

		setMultithreading(settings.numThreads);

		vlmData.load("data/MediumData.json");
	}

	void setMultithreading(unsigned int _numThreads)
	{
		// calculate number of threads according to available processors
		numThreads = max(1, min(_numThreads, numProcs));

		// create threads and samplers for each thread
		threads = new std::thread * [numThreads];
		samplers = new MTRandom * [numThreads];

		// assign different seeds to each sampler
		// Linear Congruential Generator used for seed
		// x + 1 = [a * (x - 1) + c] % m
		// where a = 48271, c = 0
		int m = pow(2, 32) - 1;
		for (unsigned int i = 0; i < numThreads; i++)
			samplers[i] = new MTRandom((48271 * (i + 1)) % m);

		// calculate number of tiles
		totalXTiles = (canvas->getWidth() + tileSize - 1) / tileSize;
		float totalYTiles = (canvas->getHeight() + tileSize - 1) / tileSize;
		totalTiles = totalXTiles * totalYTiles;

		// create samples for each tile
		tileSamples.resize(totalTiles, 1);

		tileCounter.store(0);
	}

	void clear()
	{
		film->clear();
	}

	// power heuristic for multiple importance sampling
	float powerHeuristic(const float& pdf1, const float& pdf2)
	{
		if (pdf1 + pdf2 <= EPSILON)
			return EPSILON;

		float weight1 = pdf1 * pdf1;
		float weight2 = pdf2 * pdf2;

		float sum = weight1 + weight2;

		return weight1 / sum;
	}

	// ###############################################################################################################


	// VOLUME RENDERING ##############################################################################################

	// calculate transmittance using Beer-Lambert law
	inline Color computeTransmit(const Color& sigma_t, const float& dist) const
	{
		if (sigma_t.Lum() < EPSILON)
			return Color(1.0f);

		return (sigma_t * -dist).exp();
	}

	Color endMediumTrace(PathData& path, const ShadingData& shadingData, const Color& sigma_t)
	{
		// update path throughput
		path.pathThroughput *= computeTransmit(sigma_t, shadingData.t);

		if (shadingData.bsdf->isLight())
			return misLight(shadingData, path);

		// compute surface interaction
		bool reflect = calculateSurfaceInteraction(path, shadingData);

		// check if ray is reflected, if reflected trace medium else trace surface
		if (reflect)
		{
			Color direct = computeDirect(shadingData, path.sampler, true) * path.pathThroughput;
			return direct + pathTraceMedium(path, shadingData);
		}
		else
		{
			path.depth++;					// increase path depth
			return pathTrace(path);
		}
	}

	Color computeInScattering(const Vec3& x, const PathData& path, int nLightSamples)
	{
		Color Li(0.0f);

		for (int s = 0; s < nLightSamples; s++)
		{
			LightSample light = scene->sampleLightPoint(path.sampler);
			if (light.isNull) continue;

			float lightPdf = light.pdf * light.pmf;

			Vec3 lightDir = light.p - x;
			float lightDist = lightDir.length();
			Vec3 wi = lightDir / lightDist;

			// Visibility check
			Ray emittedRay(x, wi);
			IntersectionData intersection = scene->traverse(emittedRay);

			float vlmDist = lightDist;
			float F = 0.0f;

			if (intersection.t < lightDist - EPSILON)
			{
				ShadingData shadingData = scene->calculateShadingData(intersection, emittedRay);

				if (!(shadingData.bsdf->isVolume() && scene->visible(shadingData.x, light.p)))
					continue;

				// calculate fresnel if intersected surface is dielectric
				F = ShadingHelper::fresnelDielectric(fabsf(shadingData.sNormal.dot(wi)), 1.f / path.eta);

				vlmDist = intersection.t;
			}

			// Geometry / light contribution
			float gTerm = light.isArea ? max(-Dot(lightDir, light.n), 0.0f) / (lightDist * lightDist) : 1.0f;
			Color surfaceContrib = light.emitted * gTerm * (1.f - F) / lightPdf;

			// Transmittance
			Color transmittance = computeTransmit(path.sigma_t, vlmDist);

			// Phase function
			float phaseVal = path.scatterPhaseFunc->evaluate(-path.r.dir, wi);

			// MIS weight
			float misWeight = powerHeuristic(lightPdf * nLightSamples, phaseVal);

			// Add to accumulator
			Li += transmittance * surfaceContrib * phaseVal * misWeight;
		}

		// Average
		return Li / float(nLightSamples);
	}


	Color traceMedium(PathData& path)
	{
		// get intersection data
		IntersectionData intersection = scene->traverse(path.r);

		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, path.r);

			// Analytic free-path sampling for step (from exponential distribution)
			float step = min(-logf(path.sampler->next()) / path.sigma_t.Max(), shadingData.t);

			if (step > shadingData.t - EPSILON)
				return endMediumTrace(path, shadingData, path.sigma_t);

			// russian roulette
			float russianRouletteProbability = min(path.pathThroughput.Lum(), 0.9f);
			if (russianRouletteProbability < path.sampler->next())
				return endMediumTrace(path, shadingData, path.sigma_t);
			path.pathThroughput /= russianRouletteProbability;

			// sample new scattering event
			vlmData.sample(path.sampler->next(), path.sigma_s, path.scatterPhaseFunc);

			// update path throughput with transmittance
			path.pathThroughput *= computeTransmit(path.sigma_t, step);

			if (path.pathThroughput.Lum() < EPSILON)
				return Color(0.0f);

			Vec3 x = path.r.at(step);

			// calculate in-scattering from lights
			Color inScattered = computeInScattering(x, path, 5) * path.pathThroughput;

			// Sample phase function
			Vec3 wi = path.r.dir;
			float phaseVal = path.scatterPhaseFunc->sample(path.r.dir, wi, path.sampler);
			float phasePdf = path.scatterPhaseFunc->evaluate(path.r.dir, wi);
			float phase = phaseVal / phasePdf;		// for HG and vMF the phasePdf is same as phaseVal so the phase = 1

			// update previous specular and pdf for mis
			path.prevPdf = phasePdf;
			path.prevSpecular = false;

			// calculate scattered color
			Color mediumCol = path.sigma_s / path.sigma_t;	// albedo of medium
			Color scattered = mediumCol * phase;

			// update path throughput with scattered color
			path.pathThroughput *= scattered;

			// terminate if path throughput is too low
			if (path.pathThroughput.Lum() < EPSILON)
				return inScattered;

			// Update ray for the next bounce
			path.r.init(x, wi);

			return inScattered + traceMedium(path);
		}

		return Color(0.0f);
	}

	// regular tracking volumetric path tracing
	Color pathTraceMedium(PathData& path, const ShadingData& shadingData, bool startInMedium = false)
	{
		// set extinction coefficient of the medium
		path.sigma_t = vlmData.getTotalExtinction();
		path.eta = vlmData.getRefractiveIndex();

		if (shadingData.bsdf->isLight())
		{
			path.pathThroughput *= computeTransmit(path.sigma_t, shadingData.t);
			return misLight(shadingData, path);
		}

		// check if entering volume
		bool isEntering = shadingData.n.dot(path.r.dir) < 0.0f && !startInMedium;

		// if entering sampe fresnel and reflect or transmit
		if (isEntering)
		{
			if (calculateSurfaceInteraction(path, shadingData))
			{
				path.depth++;				// increase path depth
				return pathTrace(path);
			}
		}

		return traceMedium(path);
	}

	//#################################################################################################################


	// PATH TRACE #####################################################################################################

	// return true if volume and transmittance is updated else false
	bool computeTransmitForDirectLighting(const Vec3& p1, const Vec3& p2, Color& transmittance, bool inVolume)
	{
		// compute direction and length
		Vec3 wi = (p2 - p1);
		float length = wi.length();
		wi = wi / length;

		//shadow ray to check volume between light and shading point
		Ray shadowRay(p1 + wi * EPSILON, wi);
		IntersectionData intersection = scene->traverse(shadowRay);
		if (intersection.t < length - EPSILON)
		{
			ShadingData shadingData2 = scene->calculateShadingData(intersection, shadowRay);
			if (shadingData2.bsdf->isVolume())
			{
				// update transmittance if in volume else set inVolume to true
				if (inVolume)
				{
					transmittance *= computeTransmit(vlmData.getTotalExtinction(), shadingData2.t);
					if (transmittance.Lum() < EPSILON)
						return false;
					inVolume = false;
				}
				else
					inVolume = true;

				return computeTransmitForDirectLighting(shadingData2.x, p2, transmittance, inVolume);
			}
			else
				return false;
		}

		return true;
	}

	Color computeDirect(ShadingData shadingData, Sampler* sampler, bool inVolume = false)
	{
		if (shadingData.bsdf->isPureSpecular())
		{
			return Color(0.0f, 0.0f, 0.0f);
		}

		// Light sampling part
		LightSample light = scene->sampleLightPoint(sampler);
		if (light.isNull)
			return Color(0.0f, 0.0f, 0.0f);

		Vec3 wi = light.p - shadingData.x;
		float lengthSq = light.isArea ? wi.lengthSq() : 1.0f;
		wi = wi.normalize();

		float cosThetaSurface = max(Dot(wi, shadingData.sNormal), 0.0f);
		float gTerm = cosThetaSurface / lengthSq;

		if (light.isArea)
		{
			float cosThetaLight = max(-Dot(wi, light.n), 0.0f);
			gTerm *= cosThetaLight;
		}

		if (gTerm > 0.0f)
		{
			Color transmittance(1.0f);
			// if not visible check for volume in between
			if (!scene->visible(shadingData.x, light.p))
				if (!computeTransmitForDirectLighting(shadingData.x, light.p, transmittance, inVolume))
					return Color(0.0f);

			// evaluate bsdf
			Color bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
			float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);

			// mis weight
			float lightPdf = light.pdf * light.pmf;
			float misWeight = powerHeuristic(lightPdf, bsdfPdf);

			return (transmittance * light.emitted * bsdfVal * gTerm * misWeight) / lightPdf;
		}

		return Color(0.0f);
	}

	// returns true if ray is reflected else returns false
	bool calculateSurfaceInteraction(PathData& path, const ShadingData& shadingData)
	{
		// sample surface bsdf
		Color bsdfVal;
		float bsdfPdf;
		Vec3 wi = shadingData.bsdf->sample(shadingData, path.sampler, bsdfVal, bsdfPdf);

		path.pathThroughput *= bsdfVal *
			fabsf(wi.dot(shadingData.sNormal)) / bsdfPdf;			// update path throughput

		path.prevSpecular = shadingData.bsdf->isPureSpecular();		// update previous specular flag
		path.prevPdf = shadingData.bsdf->PDF(shadingData, wi);	// update previous bsdf pdf

		// update ray direction
		path.r.init(shadingData.x + wi * EPSILON, wi);

		// check if ray is reflected
		return shadingData.sNormal.dot(wi) > 0;
	}

	Color misLight(const ShadingData& shadingData, const PathData& path)
	{
		Color emitted = shadingData.bsdf->emit(shadingData, path.r.dir);
		if (path.prevSpecular)
			return path.pathThroughput * emitted;

		float lightPdf = scene->getLightPdf(shadingData.lightIndex, -path.r.dir);
		float weight = powerHeuristic(path.prevPdf, lightPdf);

		return path.pathThroughput * emitted * weight;
	}

	Color pathTrace(PathData& path, bool startInMedium = false)
	{
		IntersectionData intersection = scene->traverse(path.r);
		ShadingData shadingData = scene->calculateShadingData(intersection, path.r);

		if (shadingData.t < FLT_MAX)
		{
			if (startInMedium || shadingData.bsdf->isVolume())
			{
				if (path.depth < settings.maxBounces)
					return pathTraceMedium(path, shadingData, startInMedium);
			}

			if (shadingData.bsdf->isLight())
				return misLight(shadingData, path);

			// calculate direct lighting
			Color direct = path.pathThroughput * computeDirect(shadingData, path.sampler);

			// max depth reached
			if (path.depth >= settings.maxBounces)
				return direct;

			// russian roulette
			float russianRouletteProbability = min(path.pathThroughput.Lum(), 0.9f);
			if (path.sampler->next() >= russianRouletteProbability)
				return direct;
			path.pathThroughput /= russianRouletteProbability;	// update path throughput

			if (path.pathThroughput.Lum() < EPSILON)
				return direct;	// terminate path if throughput is too low

			calculateSurfaceInteraction(path, shadingData);				// sample surface bsdf
			path.depth++;												// increase path depth

			return direct + pathTrace(path, false);
		}

		Color emitted = scene->background->evaluate(path.r.dir);

		if (path.prevSpecular)
			return path.pathThroughput * emitted;

		float lightPdf = scene->getLightPdf(path.r.dir);
		float misWeight = powerHeuristic(path.prevPdf, lightPdf);

		return path.pathThroughput * emitted * misWeight;
	}

	Color pathTrace(Ray& r, Sampler* sampler)
	{
		bool startInMedium = scene->inVolume(r);
		PathData pathData(r, sampler);
		return pathTrace(pathData, startInMedium);
	}

	// ###################################################################################################################

	Color direct(Ray& r, Sampler* sampler)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			// sample light
			bool inVolume = scene->inVolume(r);
			return computeDirect(shadingData, sampler, inVolume);
		}
		return scene->background->evaluate(r.dir);
	}

	Color albedo(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
				return shadingData.bsdf->emit(shadingData, shadingData.wo);

			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(r.dir);
	}

	Color viewNormals(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);
			// if (shadingData.n.dot(r.dir) < 0.0f)
			return Color(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return scene->background->evaluate(r.dir);
	}

	void lightDebug(Sampler* sampler)
	{
		if (settings.debugVMF)
			return;

		LightSample light = scene->sampleLightPoint(sampler);

		if (!light.isNull)
			drawPoint(light.p, light.emitted / (light.pdf * light.pmf));
	}

	Color viewTriangle(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			// triangle outline
			if (scene->triangles[intersection.ID].isOnEdge(r.at(intersection.t), settings.lineWidth * intersection.t))
				return Color(0.f, 1.f, 0.f);

			// albedo
			ShadingData shadingData = scene->calculateShadingData(intersection, r);

			if (shadingData.bsdf->isLight())
				return shadingData.bsdf->emit(shadingData, shadingData.wo);

			//if (shadingData.n.dot(r.dir) > 0)
			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(r.dir);
	}


	// PHASE SAMPLING TEST ##################################################################################

	void phaseTest(Sampler* sampler)
	{
		if (!settings.debugVMF)
			return;

		const float sphereRadius = use<SceneBounds>().sceneRadius;	// Radius of the sphere

		Color sigma_s;
		std::unique_ptr<PhaseFunction> phaseF;
		float pdf = vlmData.sample(sampler->next(), sigma_s, phaseF);

		Vec3 sample;
		float phasePdf = phaseF->sample(V3_UP, sample, sampler);
		float phaseVal = phaseF->evaluate(V3_UP, sample);

		std::string name = typeid(*phaseF).name();
		std::string typeidName = typeid(HGPhase).name();
		Color col = name == typeidName ? COLOR_GREEN : COLOR_RED;
		col = col.blend(name == typeidName ? COLOR_YELLOW : COLOR_BLUE, pdf * phasePdf);

		float length = sphereRadius * phasePdf * pdf;

		sample = sample * length + use<SceneBounds>().sceneCentre; // Translate sample to scene center
		drawPoint(sample, col);
	}

	// ###################################################################################################################

	void drawPoint(Vec3 point, Color col = COLOR_RED)
	{
		Vec3 dir = (point - scene->camera.origin).normalize();
		bool isInFront = dir.dot(scene->camera.viewDirection) > 0.0f;

		if (isInFront)
		{
			float px, py;
			if (scene->camera.projectOntoCamera(point, px, py))
				film->splat(px, py, col);
		}
	}

	void render()
	{
		film->incrementSPP();

		for (unsigned int y = 0; y < film->height; y++)
		{
			for (unsigned int x = 0; x < film->width; x++)
			{
				if (!settings.render)
				{
					float px = x + 0.5f;
					float py = y + 0.5f;
					Ray ray = scene->camera.generateRay(px, py);

					Color col;
					switch (settings.drawMode)
					{
					case DM_ALBEDO:
						col = albedo(ray);
						break;
					case DM_NORMALS:
						col = viewNormals(ray);
						break;
					case DM_DIRECT:
						col = direct(ray, samplers[0]);
						break;
					case DM_TRIANGLES:
						col = viewTriangle(ray);
						break;
					}

					film->splat(px, py, col);
				}
				else
				{

					float px = x + samplers[0]->next();
					float py = y + samplers[0]->next();
					Ray ray = scene->camera.generateRay(px, py);

					Color col = pathTrace(ray, samplers[0]);
					film->splat(px, py, col);
				}
			}
		}
	}

	// TILE BASED ADAPTIVE SAMPLING ########################################################################################

	void calculateTileSamples()
	{
		for (unsigned int i = 0; i < totalTiles; i++)
		{
			unsigned int startx = (i % totalXTiles) * tileSize;
			unsigned int starty = (i / totalXTiles) * tileSize;

			unsigned int endx = min(startx + tileSize, film->width);
			unsigned int endy = min(starty + tileSize, film->height);
			std::vector<float> lums = film->getLums(startx, starty, endx, endy);

			// Compute average luminance
			float total = 0.0f;
			for (float lum : lums)
				total += lum;

			float mean = (lums.empty()) ? 0.0f : total / lums.size();

			// Compute variance
			float variance = 0.0f;
			for (float lum : lums)
				variance += (lum - mean) * (lum - mean);

			variance = (lums.empty()) ? 0.0f : variance / lums.size();
			float weight = clamp(variance / (variance + mean * mean + EPSILON), EPSILON, 1.0f); // Example weighting formula

			tileSamples[i] = (unsigned int)(film->SPP + (settings.totalSPP - film->SPP) * weight);
		}
	}

	// #####################################################################################################################

	// MULTI THREADING #####################################################################################################

	void processTile(unsigned int id)
	{
		unsigned int i;
		while ((i = tileCounter.fetch_add(1)) < totalTiles)
		{
			// skip tile if it has enough samples already
			if (film->SPP > settings.initSPP && film->SPP > tileSamples[i])
				continue;

			unsigned int startx = (i % totalXTiles) * tileSize;
			unsigned int starty = (i / totalXTiles) * tileSize;

			unsigned int endx = min(startx + tileSize, film->width);
			unsigned int endy = min(starty + tileSize, film->height);

			for (unsigned int y = starty; y < endy; y++)
			{
				for (unsigned int x = startx; x < endx; x++)
				{
					if (settings.drawMode == DM_LIGHTS)
					{
						phaseTest(samplers[id]);
						lightDebug(samplers[id]); // draw light for each thread
						continue;
					}

					if (!settings.render)
					{
						float px = x + 0.5f;
						float py = y + 0.5f;
						Ray ray = scene->camera.generateRay(px, py);

						Color col;
						switch (settings.drawMode)
						{
						case DM_ALBEDO:
							col = albedo(ray);
							break;
						case DM_NORMALS:
							col = viewNormals(ray);
							break;
						case DM_DIRECT:
							col = direct(ray, samplers[id]);
							break;
						case DM_TRIANGLES:
							col = viewTriangle(ray);
							break;
						}

						film->splat(px, py, col);
					}
					else
					{
						float px = x + samplers[id]->next();
						float py = y + samplers[id]->next();
						Ray ray = scene->camera.generateRay(px, py);

						Color col = pathTrace(ray, samplers[id]);
						film->splat(px, py, col);
					}
				}
			}
		}
	}

	void renderMT()
	{
		if (settings.debug || !settings.render)
			film->clear();
		else
			film->incrementSPP();

		tileCounter.store(0);

		// process all tiles
		for (int i = 0; i < numThreads; i++)
			threads[i] = new std::thread(&Renderer::processTile, this, i);

		for (int i = 0; i < numThreads; i++)
		{
			threads[i]->join();
			delete threads[i];
		}

		// calculate samples of tiles for adaptive sampling
		if (film->SPP == settings.initSPP)
			calculateTileSamples();
	}

	// ##################################################################################################################

	void createAOV(AOV& aov)
	{
		aov = AOV(film->width, film->height);

		int sppY, spp, sppIndex;
		for (unsigned int y = 0; y < aov.height; y++)
		{
			sppY = (y / tileSize) * totalXTiles;
			for (unsigned int x = 0; x < aov.width; x++)
			{
				// calculate index
				unsigned int index = y * aov.width + x;

				sppIndex = sppY + x / tileSize;
				spp = settings.adaptiveSampling && settings.initSPP < film->SPP ? min(tileSamples[sppIndex], film->SPP) : film->SPP;

				// set colour
				Color col = film->film[index] / (float)spp;
				memcpy(&aov.color[index * 3], &col.rgb, sizeof(float) * 3);

				// create ray
				float px = x + 0.5f;
				float py = y + 0.5f;
				Ray ray = scene->camera.generateRay(px, py);

				// set albedo 
				col = albedo(ray);
				memcpy(&aov.albedo[index * 3], &col.rgb, sizeof(float) * 3);

				// set normals
				col = viewNormals(ray);
				memcpy(&aov.normal[index * 3], &col.rgb, sizeof(float) * 3);
			}
		}
	}

	void draw()
	{
		unsigned char r, g, b;
		int sppY, spp, sppIndex;
		for (unsigned int y = 0; y < film->height; y++)
		{
			sppY = (y / tileSize) * totalXTiles;
			for (unsigned int x = 0; x < film->width; x++)
			{
				sppIndex = sppY + x / tileSize;
				spp = settings.adaptiveSampling && settings.initSPP < film->SPP ? min(tileSamples[sppIndex], film->SPP) : film->SPP;
				spp = clamp(spp, 1, spp);

				film->tonemap(x, y, r, g, b, spp, settings.toneMap);
				canvas->draw(y * film->width + x, r, g, b);
			}
		}
	}

	void draw(const AOV& aov)
	{
		Color col;
		unsigned char r, g, b;
		unsigned int index, total = film->height * film->width;

		for (unsigned int i = 0; i < total; i++)
		{
			index = i * 3;
			film->tonemap(aov.output[index],
				aov.output[index + 1],
				aov.output[index + 2],
				r, g, b, settings.toneMap);
			canvas->draw(i, r, g, b);
		}
	}

	int getSPP() const
	{
		return film->SPP;
	}

	void saveHDR(std::string filename)
	{
		film->save(filename);
	}

	void savePNG(std::string filename)
	{
		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
	}

	void toggleRender()
	{
		clear();
		settings.render = !settings.render;
	}

	void toggleLightVMF()
	{
		if (settings.drawMode == DM_LIGHTS)
			settings.debugVMF = !settings.debugVMF;
	}

	void cycleDrawMode()
	{
		if (settings.render)
			return;

		clear();

		// cycle through the draw modes
		switch (settings.drawMode)
		{
		case DM_NORMALS:
			settings.drawMode = DM_TRIANGLES;
			break;
		case DM_TRIANGLES:
			settings.drawMode = DM_ALBEDO;
			break;
		case DM_ALBEDO:
			settings.drawMode = DM_LIGHTS;
			break;
		case DM_LIGHTS:
			settings.drawMode = DM_DIRECT;
			break;
		default:
			settings.drawMode = DM_NORMALS;
			break;
		}
	}
};