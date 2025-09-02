#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Bvh.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Camera.h"

class Scene
{
public:
	std::vector<Triangle> triangles;
	std::vector<BSDF*> materials;
	std::vector<Light*> lights;
	Light* background = NULL;
	BVHTree bvh;

	SceneCamera cameraMobile;
	SceneCamera cameraFixed;

	SceneCamera* camera;

	AABB bounds;

	Scene()
	{
		camera = &cameraMobile;
	}

	~Scene()
	{
		delete background;

		for (BSDF* mat : materials)
			if (mat != nullptr)
				delete mat;
	}

	void build()
	{
		// Add BVH building code here
		bvh.build(triangles, bounds);

		// Do not touch the code below this line!
		// Build light list
		for (int i = 0; i < triangles.size(); i++)
		{
			if (materials[triangles[i].materialIndex]->isLight())
			{
				AreaLight* light = new AreaLight();
				light->triangle = &triangles[i];
				light->emission = materials[triangles[i].materialIndex]->emission;
				lights.push_back(light);
			}
		}
	}

	void switchCamera(bool mobile)
	{
		if (mobile) camera = &cameraMobile;
		else camera = &cameraFixed;
	}

	IntersectionData traverseAll(const Ray& ray)
	{
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		for (int i = 0; i < triangles.size(); i++)
		{
			float t;
			float u;
			float v;
			if (triangles[i].rayIntersect(ray, t, u, v))
			{
				if (t < intersection.t)
				{
					intersection.t = t;
					intersection.ID = i;
					intersection.alpha = u;
					intersection.beta = v;
					intersection.gamma = 1.0f - (u + v);
				}
			}
		}
		return intersection;
	}

	IntersectionData traverse(const Ray& ray)
	{
		//return traverseAll(ray);
		return bvh.traverse(ray);
	}

	Light* sampleLight(Sampler* sampler, float& pmf) const
	{
		if (lights.size() <= 0)
		{
			pmf = 1.0f;		// no lights, return uniform distribution
			return nullptr;
		}

		pmf = 1 / (float)lights.size();		// probability mass function

		int i = sampler->next() * lights.size(); // sample a light index uniformly
		if (i >= lights.size())
			i = lights.size() - 1;	// ensure index is within bounds

		return lights[i];
	}

	LightSample sampleLightPoint(Sampler* sampler) const
	{
		LightSample sample;
		Light* light = sampleLight(sampler, sample.pmf);

		sample.isNull = light == nullptr;
		if (!sample.isNull)
		{
			sample.isArea = light->isArea();
			sample.p = light->sample(sampler, sample.emitted, sample.pdf);

			// adjust position if not an area light (i.e background light)
			if (!sample.isArea)
				sample.p = use<SceneBounds>().sceneCentre + (sample.p * use<SceneBounds>().sceneRadius);

			sample.n = light->normal(sample.p);
		}

		return sample;
	}

	// Do not modify any code below this line
	void init(std::vector<Triangle> meshTriangles, std::vector<BSDF*> meshMaterials, Light* _background)
	{
		for (int i = 0; i < meshTriangles.size(); i++)
		{
			triangles.push_back(meshTriangles[i]);
			bounds.extend(meshTriangles[i].vertices[0].p);
			bounds.extend(meshTriangles[i].vertices[1].p);
			bounds.extend(meshTriangles[i].vertices[2].p);
		}
		for (int i = 0; i < meshMaterials.size(); i++)
		{
			materials.push_back(meshMaterials[i]);
		}
		background = _background;
		if (background->totalIntegratedPower() > 0)
		{
			lights.push_back(background);
		}
	}
	bool visible(const Vec3& p1, const Vec3& p2)
	{
		Ray ray;
		Vec3 dir = p2 - p1;
		float maxT = dir.length() - (2.0f * EPSILON);
		dir = dir.normalize();
		ray.init(p1 + (dir * EPSILON), dir);
		return bvh.traverseVisible(ray, maxT);
	}

	Color emit(Triangle* light, ShadingData shadingData, Vec3 wi)
	{
		return materials[light->materialIndex]->emit(shadingData, wi);
	}

	ShadingData calculateShadingData(IntersectionData intersection, Ray& ray)
	{
		ShadingData shadingData = {};
		if (intersection.t < FLT_MAX)
		{
			shadingData.x = ray.at(intersection.t);
			shadingData.n = triangles[intersection.ID].n;
			shadingData.gNormal = triangles[intersection.ID].gNormal();
			triangles[intersection.ID].interpolateAttributes(intersection.alpha, intersection.beta, intersection.gamma, shadingData.sNormal, shadingData.tu, shadingData.tv);
			shadingData.bsdf = materials[triangles[intersection.ID].materialIndex];
			shadingData.wo = -ray.dir;
			if (shadingData.bsdf->isTwoSided())
			{
				if (Dot(shadingData.wo, shadingData.sNormal) < 0)
				{
					shadingData.sNormal = -shadingData.sNormal;
				}
				if (Dot(shadingData.wo, shadingData.gNormal) < 0)
				{
					shadingData.gNormal = -shadingData.gNormal;
				}
			}
			shadingData.frame.fromVector(shadingData.sNormal);
			shadingData.t = intersection.t;
		}
		else
		{
			shadingData.wo = -ray.dir;
			shadingData.t = intersection.t;
		}
		return shadingData;
	}
};