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
	std::vector<AABB> volumes;
	std::vector<BSDF*> materials;
	std::vector<Light*> lights;
	Light* background = NULL;
	BVHTree bvh;
	SceneCamera camera;
	AABB bounds;

	~Scene()
	{
		delete background;

		for (BSDF* mat : materials)
			if (mat != nullptr)
				delete mat;
	}

	void addVolume(const std::vector<Triangle>& triangles)
	{
		AABB volumeBounds;
		for (const Triangle& t : triangles)
		{
			volumeBounds.extend(t.vertices[0].p);
			volumeBounds.extend(t.vertices[1].p);
			volumeBounds.extend(t.vertices[2].p);
			this->triangles.push_back(t);
		}
		volumes.push_back(volumeBounds);
	}

	bool inVolume(const Ray& ray)
	{
		for (const AABB& volume : volumes)
			if (volume.contains(ray.o))
				return true;
		return false;
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
				light->triangle->lightIndex = lights.size(); // Set the light index for the triangle
				light->emission = materials[triangles[i].materialIndex]->emission;
				lights.push_back(light);
			}
		}
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
		// return traverseAll(ray);
		return bvh.traverse(ray);
	}

	Light* sampleLight(Sampler* sampler, float& pmf) const
	{
		pmf = 1 / (float)lights.size();		// probability mass function
		// unsigned int i = min(sampler->next() * lights.size(), lights.size() - 1.0f);
		int i = sampler->next(lights.size());
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

			// get the normal at the sample position
			// if area light , use the triangle normal
			sample.n = light->normal(sample.p.normalize());
		}

		return sample;
	}

	float getLightPdf(const int& lightIndex, const Vec3& wi) const
	{
		return (lightIndex >= 0 && lightIndex < lights.size()) ?
			lights[lightIndex]->PDF(wi) / lights.size() : 0.0f;
	}

	float getLightPdf(const Vec3& wi) const {
		return background ? background->PDF(wi) / lights.size() : 0.0f;
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
			shadingData.lightIndex = triangles[intersection.ID].lightIndex;
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