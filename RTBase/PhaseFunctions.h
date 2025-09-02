#pragma once

#include "Sampling.h"
#include "Core.h"

// phase function interface
class PhaseFunction
{
public:
	PhaseFunction() = default;
	virtual float sample(const Vec3& wo, Vec3& wi, Sampler* sampler) const = 0;
	virtual float evaluate(const Vec3& wo, const Vec3& wi) const = 0;
	virtual PhaseFunction* clone() const = 0;
	virtual void print() const = 0;
};

// Henyey-Greenstein phase function implementation
class HGPhase : public PhaseFunction
{
	float g; // Anisotropy parameter
public:
	HGPhase(float g = 1.0f) : g(g) {}

	void SetG(float _g) { g = _g; }

	float evaluate(const Vec3& wo, const Vec3& wi) const override
	{
		float cosTheta = wo.dot(wi);
		float denom = 1.0f + g * g - 2.0f * g * cosTheta;
		if (denom <= 0.0f) {
			// Should not happen with g in [-1, 1], but as a safeguard.
			return 0.0f;
		}
		return (1.0f - g * g) / (4.0f * M_PI * denom * sqrtf(denom));
	}

	float sample(const Vec3& wo, Vec3& wi, Sampler* sampler) const override
	{
		float r1, r2;
		r1 = sampler->next();
		r2 = sampler->next();

		float cosTheta;
		if (fabsf(g) < 1e-3) {
			// For very small g, the distribution is nearly isotropic.
			// Sample a uniform sphere to avoid numerical issues.
			cosTheta = 1.0f - 2.0f * r1;
		}
		else {
			// Standard HG sampling formula
			float g2 = g * g;
			float term = (1.0f - g2) / (1.0f - g + 2.0f * g * r1);
			cosTheta = (1.0f + g2 - term * term) / (2.0f * g);
		}

		float sinTheta = sqrtf(max(0.0f, 1.0f - cosTheta * cosTheta));
		float phi = 2.0f * M_PI * r2;

		// Convert spherical coordinates to a vector in a local coordinate system
		Vec3 localDir(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);

		// create a frame based on the outgoing direction 'wo'
		Frame frame;
		frame.fromVector(wo);
		wi = frame.toWorld(localDir);

		// The PDF for the sampled direction is simply the evaluation of the phase function.
		return evaluate(wo, wi);
	}

	void print() const override {
		std::cout << "HG with g = " << g << std::endl;
	}

	HGPhase* clone() const override {
		return new HGPhase(g);
	}

};

// von Mises-Fisher distribution data structure
struct vMFData
{
	Vec3 mu;			// Mean direction (unit vector)
	float kappa;		// Concentration parameter
	float weight;		// Weight of the sample

	// Precomputed values
	float denom;
	float norm;
};

// von Mises-Fisher phase function implementation
class vMFPhase : public PhaseFunction
{
	std::vector<vMFData> samples;	// Loaded vMF components
	std::vector<float> cumPDF;		// Cumulative weights

	float calculateIndividualPdf(const vMFData& data, const Vec3& dir) const
	{
		if (data.kappa < 1e-4f)
			return 1.0f / (4.0f * float(M_PI));

		float dot = data.mu.dot(dir);
		return data.norm * std::exp(data.kappa * dot);
	}

	float calculatePdf(const Vec3& wi) const
	{
		if (samples.empty())
			return 0.0f;

		float pdf = 0.0f;
		for (const auto& comp : samples)
			pdf += comp.weight * calculateIndividualPdf(comp, wi);
		return pdf;
	}

	void calculateCumulativePDF()
	{
		cumPDF.resize(samples.size());
		float sum = 0.0f;
		for (size_t i = 0; i < samples.size(); ++i)
		{
			sum += samples[i].weight;
			cumPDF[i] = sum;
		}
		if (sum > 0.0f)
		{
			for (auto& pdf_val : cumPDF)
				pdf_val /= sum;

			for (auto& data : samples)
				data.weight /= sum;
		}
	}

	Vec3 sampleDirection(float kappa, const Vec3& mu, float r1, float r2) const
	{
		float cosTheta;
		if (kappa < 1e-5f)
		{
			cosTheta = 1.0f - 2.0f * r1;
		}
		else
		{
			float logTerm = std::log(1.0f - r1 + r1 * std::exp(-2.0f * kappa));
			cosTheta = 1.0f + (1.0f / kappa) * logTerm;
		}
		float sinTheta = std::sqrt(fmax(0.0f, 1.0f - cosTheta * cosTheta));
		float phi = 2.0f * float(M_PI) * r2;

		Vec3 local(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);

		Frame frame;
		frame.fromVector(mu);
		return frame.toWorld(local);
	}
public:
	vMFPhase() = default;

	vMFPhase(const std::vector<vMFData>& samples)
	{
		this->samples = samples;

		for (auto& data : this->samples)
		{
			// Precompute the denominator for PDF calculation
			data.denom = fmax(1e-6f, 4.0f * float(M_PI) * std::sinh(data.kappa));
			data.norm = data.kappa / data.denom;
		}

		calculateCumulativePDF();
	}

	float sample(const Vec3& wo, Vec3& wi, Sampler* sampler) const override
	{
		float r1 = sampler->next();
		float r2 = sampler->next();
		float r3 = sampler->next();

		if (samples.empty())
			return 0.0f;

		size_t index = 0;
		if (cumPDF.size() > 1)
		{
			auto it = std::lower_bound(cumPDF.begin(), cumPDF.end(), r1);
			index = std::distance(cumPDF.begin(), it);
			if (index >= samples.size()) index = samples.size() - 1;
		}

		const vMFData& data = samples[index];
		Vec3 localWi = sampleDirection(data.kappa, data.mu, r2, r3);

		Frame frame;
		frame.fromVector(wo);
		wi = frame.toWorld(localWi); // Convert to world coordinates

		return calculatePdf(localWi);
	}

	float evaluate(const Vec3& wo, const Vec3& wi) const override
	{
		if (samples.empty())
			return 0.0f;

		Frame frame;
		frame.fromVector(wo); // Create frame from outgoing direction
		Vec3 localWi = frame.toLocal(wi); // Convert incoming direction to local coordinates

		return calculatePdf(localWi);
	}

	void print() const override {
		std::cout << "vMFDistribution with vmfs \n";
		for (const auto& vmf : samples)
			std::cout << "  mu = " << vmf.mu << ", kappa = " << vmf.kappa << ", weight = " << vmf.weight << "\n";
	}

	vMFPhase* clone() const override {
		return new vMFPhase(samples);
	}
};