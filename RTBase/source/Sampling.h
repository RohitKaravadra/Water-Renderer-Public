#pragma once

#include "Core.h"
#include <random>
#include <algorithm>

const float Inv4Pi = 0.07957747154594767f; // 1 / (4 * M_PI)

class Sampler
{
public:
	virtual float next() = 0;
	virtual int next(const int& range) = 0;
};

class MTRandom : public Sampler
{
public:
	std::mt19937 generator;
	std::uniform_real_distribution<float> dist;

	MTRandom(unsigned int seed = 1) : dist(0, 1.0f)
	{
		generator.seed(seed);
	}

	float next() override
	{
		return std::max(dist(generator), 1e-8f);
	}

	// smaple an integer in the range [min, max)
	int next(const int& range) override
	{
		return static_cast<int>(next() * range);
	}
};

// Note all of these distributions assume z-up coordinate system
class SamplingDistributions
{
public:
	static Vec3 uniformSampleHemisphere(float r1, float r2)
	{
		float cosTheta = r1;
		float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
		float phi = 2.0f * M_PI * r2;
		return Vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
	}

	static float uniformHemispherePDF(const Vec3& wi)
	{
		return fabsf(wi.z) / (2.0f * M_PI);
	}

	static Vec3 cosineSampleHemisphere(float r1, float r2)
	{
		float cosTheta = sqrtf(r1);
		float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
		float phi = 2.0f * M_PI * r2;
		return Vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
	}

	static float cosineHemispherePDF(const Vec3& wi)
	{
		return fabsf(wi.z) / M_PI;
	}

	static Vec3 uniformSampleSphere(float r1, float r2)
	{
		float cosTheta = 1.0f - 2.0f * r1;
		float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
		float phi = 2.0f * M_PI * r2;
		return Vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
	}

	static float uniformSpherePDF(const Vec3& wi)
	{
		return 1.0f / (4.0f * M_PI);
	}
};