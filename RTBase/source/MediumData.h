#pragma once

#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "json.hpp"
#include "PhaseFunctions.h"

using json = nlohmann::json;

/**
 * @brief Per-particle optical properties (absorption, scattering, extinction).
 * Scaled internally by concentration or mass density.
 */
struct ParticleProperties
{
	std::string name;
	Color sigma_a;      // effective absorption coefficient
	Color sigma_s;      // effective scattering coefficient
	Color sigma_t;      // effective extinction coefficient
	std::unique_ptr<PhaseFunction> phaseF;
	float concentration;  // number concentration (#/m^3)
	float eta;            // refractive index (real)

	ParticleProperties(const std::string& name,
		const Color& sigma_a, const Color& sigma_s,
		float concentration, float eta = 1.0f, float mass = 0.0f,
		int phaseType = -1, float g = 0.0f, const std::vector<vMFData>& vmf = {}
	)
		: name(name), sigma_a(sigma_a* concentration), sigma_s(sigma_s* concentration),
		sigma_t((sigma_a + sigma_s)* concentration), concentration(concentration),
		eta(eta)
	{
		if (phaseType == 0)			// Henyey–Greenstein
			phaseF = std::make_unique<HGPhase>(g);
		else if (phaseType == 1)	// von Mises–Fisher
			phaseF = std::make_unique<vMFPhase>(vmf);
		else						// isotropic fallback
			phaseF = std::make_unique<IsotropicPhase>();
	}
};

/**
 * @brief Represents a full participating medium (host + suspended particles).
 */
class MediumData
{
private:
	std::vector<ParticleProperties> particles;
	std::vector<float> cdf;

	// Host properties
	std::string host_name;
	Color host_sigma_a;
	Color host_sigma_s;
	float host_eta = 1.0f;

	// Combined properties
	Color sigma_t;
	float eta = 1.0f;

	void computeBulkProperties()
	{
		// Start with host medium
		sigma_t = host_sigma_a + host_sigma_s;

		// Add particle contributions
		for (const auto& p : particles)
			sigma_t += p.sigma_t;

		// --- Effective refractive index (Clausius–Mossotti approx.) ---
		// n_eff ≈ sqrt( (1 + 2f) / (1 - f) ), f = Σ φ_i (polarizability fraction)
		// Here we simplify to weighted volume fraction of particle indices
		float total_weight = 1.0f; // host counts as 1
		float weighted_eta = host_eta;

		for (const auto& p : particles)
		{
			if (p.eta > 0.0f && p.concentration > 0.0f)
			{
				weighted_eta += p.eta * p.concentration;
				total_weight += p.concentration;
			}
		}
		eta = (total_weight > 0.0f) ? weighted_eta / total_weight : host_eta;
	}

	void build()
	{
		cdf.clear();
		float total_concentration = 0.0f;

		for (const auto& d : particles)
		{
			total_concentration += d.concentration;
			cdf.emplace_back(total_concentration);
		}

		if (total_concentration > 0.0f)
		{
			for (auto& val : cdf)
				val /= total_concentration;
		}

		computeBulkProperties();
	}

	static int binarySearch(const std::vector<float>& cdf, float value)
	{
		auto it = std::lower_bound(cdf.begin(), cdf.end(), value);
		return static_cast<int>(std::distance(cdf.begin(), it));
	}

public:
	MediumData() = default;

	Color getTotalExtinction() const { return sigma_t; }
	float getRefractiveIndex() const { return eta; }

	void load(const std::string& filepath)
	{
		std::ifstream file(filepath);
		if (!file.is_open())
		{
			std::cerr << "Error: Failed to open medium data file: " << filepath << std::endl;
			return;
		}

		nlohmann::json j;
		try { file >> j; }
		catch (nlohmann::json::parse_error& e)
		{
			std::cerr << "Error: JSON parsing failed: " << e.what() << std::endl;
			return;
		}
		file.close();

		// Host
		if (j.contains("host"))
		{
			const auto& host = j["host"];
			host_name = host.value("name", "Unknown");
			host_sigma_a = { host["sigma_a"][0], host["sigma_a"][1], host["sigma_a"][2] };
			host_sigma_s = { host["sigma_s"][0], host["sigma_s"][1], host["sigma_s"][2] };
			host_eta = host.value("eta", 1.0f);
		}
		else
		{
			host_name = "Air";
			host_sigma_a = Color(0.0f);
			host_sigma_s = Color(0.0f);
			host_eta = 1.0f;
		}

		// Preload vMF sets
		std::vector<std::vector<vMFData>> vmf_definitions;
		if (j.contains("vmf"))
		{
			for (const auto& vmf_group : j["vmf"])
			{
				std::vector<vMFData> current;
				for (const auto& p : vmf_group["data"])
				{
					Vec3 mu = { p["mu"][0], p["mu"][1], p["mu"][2] };
					current.emplace_back(vMFData{ mu, p["kappa"], p["weight"] });
				}
				vmf_definitions.push_back(current);
			}
		}

		// Particles
		particles.clear();
		if (j.contains("particles"))
		{
			for (const auto& p : j["particles"])
			{
				std::string name = p.value("name", "Unnamed");
				Color sigma_a_raw = { p["sigma_a"][0], p["sigma_a"][1], p["sigma_a"][2] };
				Color sigma_s_raw = { p["sigma_s"][0], p["sigma_s"][1], p["sigma_s"][2] };
				float concentration = p.value("concentration", 0.0f);
				float mass = p.value("mass", 0.0f);
				float eta = p.value("eta", 0.0f);
				int phaseType = p.value("phase_type", -1);

				if (concentration <= 0.0f) continue;

				if (phaseType == 0) { // HG
					float g = p.value("phase", 0.0f);
					particles.emplace_back(name, sigma_a_raw, sigma_s_raw, concentration, eta, mass, 0, g);
				}
				else if (phaseType == 1) { // vMF
					int vmf_index = static_cast<int>(p.value("phase", -1.0f));
					if (vmf_index >= 0 && vmf_index < vmf_definitions.size()) {
						particles.emplace_back(name, sigma_a_raw, sigma_s_raw, concentration, eta, mass, 1, 0.0f, vmf_definitions[vmf_index]);
					}
					else {
						std::cerr << "Warning: Invalid vMF index " << vmf_index
							<< " for particle '" << name << "'. Defaulting isotropic." << std::endl;
						particles.emplace_back(name, sigma_a_raw, sigma_s_raw, concentration, eta, mass);
					}
				}
				else { // isotropic
					particles.emplace_back(name, sigma_a_raw, sigma_s_raw, concentration, eta, mass);
				}
			}
		}

		build();
		print();
	}

	void print() const
	{
		std::cout << "--- Medium Properties ---\n";
		std::cout << "Host Medium: " << host_name << "\n";
		std::cout << "  - Host Sigma_a: " << host_sigma_a << "\n";
		std::cout << "  - Host Sigma_s: " << host_sigma_s << "\n";
		std::cout << "  - Host Eta    : " << host_eta << "\n";

		std::cout << "\n--- Particles ---\n";
		for (size_t i = 0; i < particles.size(); i++)
		{
			const auto& p = particles[i];
			float pdf = (cdf.empty() || i >= cdf.size()) ? 0.0f : ((i == 0) ? cdf[0] : cdf[i] - cdf[i - 1]);

			std::cout
				<< "\n Particle        : " << p.name
				<< "\n  - Raw Sigma_a   : " << (p.sigma_a / max(p.concentration, 1e-6f))
				<< "\n  - Raw Sigma_s   : " << (p.sigma_s / max(p.concentration, 1e-6f))
				<< "\n  - Concentration : " << p.concentration
				<< "\n  - Eta           : " << p.eta
				<< "\n  - Eff Sigma_a   : " << p.sigma_a
				<< "\n  - Eff Sigma_s   : " << p.sigma_s
				<< "\n  - PDF           : " << pdf
				<< "\n  - Phase Func.   : ";
			p.phaseF->print();
		}

		std::cout << "\n\n--- Final Combined ---\n";
		std::cout << "Total Sigma_t : " << sigma_t << "\n";
		std::cout << "Final Eta     : " << eta << "\n";
		std::cout << "-------------------------\n";
	}

	float sample(float u, Color& out_sigma_s, std::unique_ptr<PhaseFunction>& out_phaseF)
	{
		if (particles.empty())
		{
			out_sigma_s = Color(0.0f);
			out_phaseF = std::make_unique<IsotropicPhase>();
			return 1.f;
		}

		int index = binarySearch(cdf, u);
		if (index < particles.size())
		{
			out_sigma_s = particles[index].sigma_s;
			out_phaseF.reset(particles[index].phaseF->clone());
		}
		else
		{
			out_sigma_s = Color(0.0f);
			out_phaseF = std::make_unique<IsotropicPhase>();
			return 1.f;
		}

		// return pdf of selecting this particle
		return index == 0 ? cdf[0] : cdf[index] - cdf[index - 1];
	}
};
