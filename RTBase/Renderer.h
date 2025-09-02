#pragma once

#include "Sampling.h"
#include "Geometry.h"
#include "Materials.h"
#include "Lights.h"
#include "GamesEngineeringBase.h"
#include "Settings.h"
#include "Debugger.h"
#include "PhaseFunctions.h"
#include "ShapeProperties.h"
#include "json.hpp"

#include <thread>
#include <functional>

using json = nlohmann::json;

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

	SETTINGS settings;						// settings for the ray tracer

	std::vector<std::pair<Vec3, float>> traceData;		// stores angle with throughputs
	std::mutex traceMutex;					// mutex for thread safety

	Debugger debugger;				// debugger for rendering debug lines

	// optical properties for heterogeneous media
	OpticalProperties optProps = Plankton_Physical;

	HGPhase ghPhase;				// Henyey-Greenstein phase function
	vMFPhase vmfFunc;				// von Mises-Fisher phase function

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
		film->init((unsigned int)scene->camera->width,
			(unsigned int)scene->camera->height, settings.filter);
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;

		setMultithreading(settings.numThreads);

		debugger.init(film, scene);

		ghPhase = HGPhase(optProps.g);
		loadVMFData("data/phase_data.json");
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

		tileCounter.store(0);
	}

	void loadVMFData(const std::string& filename)
	{
		std::vector<vMFData> data;
		std::ifstream infile(filename);

		if (!infile.is_open())
		{
			std::cerr << "Error opening file: " << filename << std::endl;
			return; // Exit if the file can't be opened
		}

		try
		{
			// 1. Parse the entire file stream into a JSON object
			json jsonData;
			infile >> jsonData;

			// 2. Check if the top-level "vmf" key exists and is an array
			if (jsonData.contains("vmf") && jsonData["vmf"].is_array())
			{
				// 3. Iterate over each element in the "vmf" array
				for (const auto& item : jsonData["vmf"])
				{
					vMFData vmf;

					// 4. Extract the mu vector, checking its format
					if (item.contains("mu") && item["mu"].is_array() && item["mu"].size() == 3)
					{
						vmf.mu.x = item["mu"][0].get<float>();
						vmf.mu.y = item["mu"][1].get<float>();
						vmf.mu.z = item["mu"][2].get<float>();
					}

					// 5. Extract kappa and weight
					if (item.contains("kappa")) {
						vmf.kappa = item["kappa"].get<float>();
					}
					if (item.contains("weight")) {
						vmf.weight = item["weight"].get<float>();
					}

					data.push_back(vmf);
				}
			}
			else
			{
				std::cerr << "Error: JSON file does not contain a 'vmf' array." << std::endl;
			}
		}
		catch (json::parse_error& e)
		{
			// Catch any errors during the parsing process
			std::cerr << "JSON parsing error: " << e.what() << std::endl;
		}

		// 6. Assign the loaded data to your global phase function object
		vmfFunc = vMFPhase(data);
	}

	void clear()
	{
		if (traceData.size() > 0)
			traceData.clear();
		film->clear();
	}

	// #################################################################################################################


	// ###################################################################################################################

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
			// return Color(shadingData.sNormal.x, shadingData.sNormal.y, shadingData.sNormal.z);
			// if (shadingData.n.dot(r.dir) < 0.0f)
			return Color(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return scene->background->evaluate(r.dir);
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

	// ###################################################################################################################


	// return true if reflected, false if transmitted
	bool sampleSurface(Ray& r, const ShadingData& shadingData, Sampler* sampler)
	{
		float pdf;
		Color bsdf;
		Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, bsdf, pdf);

		r.init(shadingData.x + (wi * EPSILON), wi);

		return  shadingData.sNormal.dot(wi) > 0;
	}

	bool findTracePoint(Ray& r, Sampler* sampler, float& pathThroughput,
		int samples = 0, bool isInShape = false, bool firstBounce = true)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);

			Vec3 prevPoint = r.o;
			bool wasInShape = isInShape;
			float prevThroughput = pathThroughput;

			if (isInShape)
			{
				// sample step from exponential distribution
				float step = -logf(1.0f - sampler->next()) / optProps.sigma_t;

				// check for exit or oversampling
				if (step > shadingData.t || step < EPSILON)
				{
					pathThroughput *= expf(-optProps.sigma_t * shadingData.t) * (optProps.sigma_s / optProps.sigma_t);
					isInShape = isInShape ^ sampleSurface(r, shadingData, sampler); // XOR
					samples = 0; // reset
				}
				else
				{
					Vec3 wi;
					float phase = ghPhase.sample(r.dir, wi, sampler);
					pathThroughput *= expf(-optProps.sigma_t * step) * (optProps.sigma_s / optProps.sigma_t);
					r.init(r.o + r.dir * step, wi);
					samples++;
				}
			}
			else
			{
				if (!sampleSurface(r, shadingData, sampler))
				{
					isInShape = true;
					samples = 0;
				}
			}

			// terminate if throughput is too low
			if (pathThroughput <= 0.0f)
				return false;

			if (!settings.render)
				debugger.addLine(DebugLine(prevPoint, r.o, (wasInShape ? COLOR_GREEN : COLOR_RED) * prevThroughput, true));

			return findTracePoint(r, sampler, pathThroughput, samples, isInShape, false);
		}

		float t;
		if (!firstBounce && !settings.render && use<SceneBounds>().bounds.rayIntersect(r, t))
			debugger.addLine(DebugLine(r.o, r.at(t), COLOR_CYAN * pathThroughput, true));

		return !firstBounce;
	}

	void shapeTrace(Sampler* sampler)
	{
		if (settings.paused)
			return;

		float sceneRadius = use<SceneBounds>().sceneRadius;
		Vec3 sceneCentre = use<SceneBounds>().sceneCentre;

		// sample a random direction from sphere
		Vec3 dir = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());

		// convert to scene point
		Vec3 p = dir * sceneRadius + sceneCentre;

		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());

		Frame frame;
		frame.fromVector(-dir);
		wi = frame.toWorld(wi);

		Ray r(p, wi);

		float pathThroughput = 1.0f;
		if (!findTracePoint(r, sampler, pathThroughput, 10, false))
			return;

		if (!settings.render)
			return;

		float t;
		if (use<SceneBounds>().bounds.rayIntersect(r, t))
		{
			Vec3 sphereP = r.at(t);

			Vec3 wo = (sphereP - sceneCentre).normalize();

			Frame frame;
			frame.fromVector(wi);
			wo = frame.toLocal(wo);

			if (settings.render)
			{
				// store trace data
				traceMutex.lock();
				traceData.push_back(std::make_pair(wo, pathThroughput));
				if (traceData.size() >= settings.maxTraces)
					saveTraceData();
				traceMutex.unlock();
			}

			frame.fromVector(V3_UP);
			wo = frame.toWorld(wo);

			sphereP = wo * use<SceneBounds>().sceneRadius;

			debugger.addLine(DebugLine(use<SceneBounds>().sceneCentre, sphereP, COLOR_CYAN * pathThroughput, true));
		}
	}

	void saveTraceData()
	{
		if (traceData.size() == 0)
		{
			std::cout << "No Trace Data recorded yet \n";
			return;
		}

		std::ofstream file("data/traces.txt", std::ios::out);
		if (!file.is_open())
		{
			std::cout << "Could not open trace file!" << std::endl;
			return;
		}

		std::cout << "Saving trace data with size " << traceData.size() << std::endl;

		for (auto& data : traceData)
			file << data.first.x << " " << data.first.y << " " << data.first.z << " " << data.second << "\n";

		file.close();

		std::cout << "Trace data saved to traces.txt" << std::endl;

		traceData.clear();
	}

	void vMFDebug(Sampler* sampler)
	{
		if (settings.paused)
			return;

		Vec3 wi;
		float phase = vmfFunc.sample(V3_UP, wi, sampler);

		float sceneRadius = use<SceneBounds>().sceneRadius;
		Vec3 sceneCentre = use<SceneBounds>().sceneCentre;

		debugger.addLine(DebugLine(sceneCentre, sceneCentre + wi * sceneRadius * phase, COLOR_WHITE, true));
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
					Ray ray = scene->camera->generateRay(px, py);

					Color col;
					switch (settings.drawMode)
					{
					case DM_ALBEDO:
						col = albedo(ray);
						break;
					case DM_NORMALS:
						col = viewNormals(ray);
						break;
					case DM_TRIANGLES:
						col = viewTriangle(ray);
						break;
					}

					film->splat(px, py, col);
				}
				else
				{

				}
			}
		}
	}

	// MULTI THREADING #####################################################################################################

	void processTile(unsigned int id)
	{
		unsigned int i;
		while ((i = tileCounter.fetch_add(1)) < totalTiles)
		{
			if (settings.render)
				continue;

			unsigned int startx = (i % totalXTiles) * tileSize;
			unsigned int starty = (i / totalXTiles) * tileSize;

			unsigned int endx = min(startx + tileSize, film->width);
			unsigned int endy = min(starty + tileSize, film->height);

			for (unsigned int y = starty; y < endy; y++)
			{
				for (unsigned int x = startx; x < endx; x++)
				{
					float px = x + 0.5f;
					float py = y + 0.5f;
					Ray ray = scene->camera->generateRay(px, py);

					Color col;
					switch (settings.drawMode)
					{
					case DM_ALBEDO:
						col = albedo(ray);
						break;
					case DM_NORMALS:
						col = viewNormals(ray);
						break;
					case DM_TRIANGLES:
						col = viewTriangle(ray);
						break;
					}

					film->splat(px, py, col);
				}
			}
		}
	}

	void renderMT()
	{
		if (!settings.render)
			film->clear();
		else
			film->incrementSPP();

		tileCounter.store(0);

		if (!settings.paused)
			debugger.clear();

		for (int i = 0; i < settings.tracePerFrame; i++)
		{
			if (settings.drawMode == DM_VMF && !settings.render)
				vMFDebug(samplers[0]);
			else
				shapeTrace(samplers[0]);
		}

		// process all tiles
		for (int i = 0; i < numThreads; i++)
			threads[i] = new std::thread(&Renderer::processTile, this, i);

		for (int i = 0; i < numThreads; i++)
		{
			threads[i]->join();
			delete threads[i];
		}

		debugger.draw();
	}

	void draw()
	{
		unsigned char r, g, b;
		int spp = clamp(film->SPP, 1, film->SPP);

		for (unsigned int y = 0; y < film->height; y++)
		{
			for (unsigned int x = 0; x < film->width; x++)
			{
				film->tonemap(x, y, r, g, b, spp, settings.toneMap);
				canvas->draw(y * film->width + x, r, g, b);
			}
		}
	}

	int getSPP() const
	{
		return film->SPP;
	}

	void toggleRender()
	{
		clear();
		settings.render = !settings.render;
		scene->switchCamera(!settings.render);
	}

	void cycleDrawMode()
	{
		clear();

		// cycle through the draw modes
		settings.drawMode = (DRAW_MODE)(((int)settings.drawMode + 1) % 4);
		scene->switchCamera(!settings.render);
	}
};