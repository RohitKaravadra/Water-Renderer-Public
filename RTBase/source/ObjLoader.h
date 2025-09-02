#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

struct OBJPoint
{
	float x, y, z;
};

struct OBJVertex {
	OBJPoint p, n; // Position and normal
};

struct OBJTriangle {
	int v1, v2, v3; // Indices of vertices
};

static class OBJLoader
{
public:

	// Load an OBJ file and populate the vertex and triangle arrays
	static bool Load(const std::string& filename, std::vector<OBJVertex>& outVertices, std::vector<OBJTriangle>& outTriangles) {
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Failed to open file: " << filename << "\n";
			return false;
		}

		std::vector<float> positions;
		std::vector<float> normals;
		std::vector<int> positionIndices, normalIndices;

		std::string line;
		while (std::getline(file, line)) {
			std::istringstream iss(line);
			std::string type;
			iss >> type;

			if (type == "v") {
				float x, y, z;
				iss >> x >> y >> z;
				positions.push_back(x);
				positions.push_back(y);
				positions.push_back(z);
			}
			else if (type == "vn") {
				float nx, ny, nz;
				iss >> nx >> ny >> nz;
				normals.push_back(nx);
				normals.push_back(ny);
				normals.push_back(nz);
			}
			else if (type == "f") {
				std::string v1, v2, v3;
				iss >> v1 >> v2 >> v3;
				std::string verts[3] = { v1, v2, v3 };

				for (std::string& v : verts) {
					size_t pos1 = v.find('/');
					size_t pos2 = v.find('/', pos1 + 1);

					int vi = std::stoi(v.substr(0, pos1)) - 1;
					int ni = std::stoi(v.substr(pos2 + 1)) - 1;

					positionIndices.push_back(vi);
					normalIndices.push_back(ni);
				}

				int idx = static_cast<int>(positionIndices.size()) - 3;
				outTriangles.emplace_back(OBJTriangle{ idx, idx + 1, idx + 2 });
			}
		}

		// Combine positions and normals into final vertex array
		for (size_t i = 0; i < positionIndices.size(); ++i) {
			int pi = positionIndices[i] * 3;
			int ni = normalIndices[i] * 3;

			OBJVertex v;
			v.p = OBJPoint{ positions[pi],positions[pi + 1],positions[pi + 2] };
			v.n = OBJPoint{ normals[ni],normals[ni + 1],normals[ni + 2] };

			outVertices.emplace_back(v);
		}

		return true;
	}
};