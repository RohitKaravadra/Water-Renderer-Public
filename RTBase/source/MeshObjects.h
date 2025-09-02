#pragma once
#include "Core.h"
#include "Geometry.h"
#include <set>

class CustomMesh
{

public:
	std::vector<Triangle> triangles;

	void invertNormals()
	{
		for (Triangle& t : triangles)
			t.invertNormals();
	}
};

class PlaneMesh : public CustomMesh
{
public:
	PlaneMesh(Vec3 center, Vec3 size, Vec3 part, unsigned int matIndex = 0)
	{
		unsigned int xPart = static_cast<unsigned int>(part.x);
		unsigned int yPart = static_cast<unsigned int>(part.y);

		Vec3 origin = center - Vec3(size.x * 0.5f, 0, -size.y * 0.5f); // Flip Z for origin
		float dx = size.x / xPart;
		float dy = size.y / yPart;

		for (unsigned int i = 0; i < xPart; i++)
		{
			for (unsigned int j = 0; j < yPart; j++)
			{
				float u0 = (float)i / xPart;
				float v0 = (float)j / yPart;
				float u1 = (float)(i + 1) / xPart;
				float v1 = (float)(j + 1) / yPart;

				// Flip Z in positions
				Vec3 p00 = origin + Vec3(i * dx, 0, -(j * dy));
				Vec3 p01 = origin + Vec3(i * dx, 0, -((j + 1) * dy));
				Vec3 p10 = origin + Vec3((i + 1) * dx, 0, -(j * dy));
				Vec3 p11 = origin + Vec3((i + 1) * dx, 0, -((j + 1) * dy));

				Vec3 normal = Vec3(0, 1, 0); // Y up, normal unchanged because it points up

				Vertex v00 = { p00, normal, u0, v0 };
				Vertex v01 = { p01, normal, u0, v1 };
				Vertex v10 = { p10, normal, u1, v0 };
				Vertex v11 = { p11, normal, u1, v1 };

				Triangle t1, t2;
				// Reverse winding order: (v0, v2, v1)
				t1.init(v00, v10, v01, matIndex);
				t2.init(v10, v11, v01, matIndex);

				triangles.push_back(t1);
				triangles.push_back(t2);
			}
		}
	}
};

class CubeMesh :public CustomMesh
{
	const int faces[6][4] = {
		{0, 1, 2, 3}, // Back (-Z)
		{5, 4, 7, 6}, // Front (+Z)
		{4, 0, 3, 7}, // Left (-X)
		{1, 5, 6, 2}, // Right (+X)
		{3, 2, 6, 7}, // Top (+Y) --> This face gets partitioned
		{4, 5, 1, 0}  // Bottom (-Y)
	};

private:

	// Simple linear interpolation between two Vec3
	Vec3 lerp(const Vec3& a, const Vec3& b, float t) const
	{
		return a + (b - a) * t;
	}

	Vec3 colToVec(const Color& col) const
	{
		return Vec3(col.r, col.g, col.b);
	}

public:
	CubeMesh(Vec3 center, Vec3 size, unsigned int matIndex = 0)
	{
		Vec3 half = size * 0.5f;

		Vec3 p[8] = {
			center + Vec3(-half.x, -half.y, -(-half.z)), // Flip Z
			center + Vec3(half.x, -half.y, -(-half.z)),
			center + Vec3(half.x, half.y, -(-half.z)),
			center + Vec3(-half.x, half.y, -(-half.z)),
			center + Vec3(-half.x, -half.y, -(half.z)),
			center + Vec3(half.x, -half.y, -(half.z)),
			center + Vec3(half.x, half.y, -(half.z)),
			center + Vec3(-half.x, half.y, -(half.z)),
		};

		Vec3 normals[6] = {
			Vec3(0, 0, 1),  // Back (-Z flipped to +Z)
			Vec3(0, 0, -1), // Front (+Z flipped to -Z)
			Vec3(-1, 0, 0), // Left (-X)
			Vec3(1, 0, 0),  // Right (+X)
			Vec3(0, 1, 0),  // Top (+Y)
			Vec3(0, -1, 0)  // Bottom (-Y)
		};

		for (int f = 0; f < 6; ++f)
		{
			Vertex v0 = { p[faces[f][0]], normals[f], 0.0f, 0.0f };
			Vertex v1 = { p[faces[f][1]], normals[f], 1.0f, 0.0f };
			Vertex v2 = { p[faces[f][2]], normals[f], 1.0f, 1.0f };
			Vertex v3 = { p[faces[f][3]], normals[f], 0.0f, 1.0f };

			Triangle t1, t2;
			// Reverse winding order
			t1.init(v0, v1, v3, matIndex);
			t2.init(v1, v2, v3, matIndex);

			triangles.push_back(t1);
			triangles.push_back(t2);
		}
	}

	// create mesh
	CubeMesh(Vec3 center, Vec3 size, Vec3 partition, unsigned int matIndex = 0, std::string removeFacesStr = "\0")
	{
		CreateMesh(center, size, partition, matIndex, "\0", "\0", 0.0f, removeFacesStr);
	}

	// create mesh with normal map
	CubeMesh(Vec3 center, Vec3 size, Vec3 partition,
		std::string normalMap, unsigned int matIndex = 0, std::string removeFacesStr = "\0")
	{
		CreateMesh(center, size, partition, matIndex, normalMap, "\0", 0.0f, removeFacesStr);
	}

	// create mesh with height map
	CubeMesh(Vec3 center, Vec3 size, Vec3 partition, std::string heightMap,
		float strength, unsigned int matIndex = 0, std::string removeFacesStr = "\0")
	{
		CreateMesh(center, size, partition, matIndex, "", heightMap, strength, removeFacesStr);
	}

	// create mesh with normal and height map
	CubeMesh(Vec3 center, Vec3 size, Vec3 partition, std::string normalMap,
		std::string heightMap, float strength, unsigned int matIndex = 0, std::string removeFacesStr = "\0")
	{
		CreateMesh(center, size, partition, matIndex, normalMap, heightMap, strength, removeFacesStr);
	}

	void CreateMesh(Vec3 center, Vec3 size, Vec3 partition, unsigned int matIndex,
		std::string normalMap, std::string heightMap, float strength,
		std::string removeFacesStr)
	{
		// --- Parse faces to remove ---
		std::set<int> removedFaces;
		if (!removeFacesStr.empty())
		{
			static std::map<std::string, int> faceMap = {
				{"back", 0}, {"front", 1}, {"left", 2},
				{"right", 3}, {"top", 4}, {"bottom", 5}
			};

			std::istringstream iss(removeFacesStr);
			std::string face;
			while (iss >> face)
			{
				std::transform(face.begin(), face.end(), face.begin(), ::tolower);
				if (faceMap.count(face))
					removedFaces.insert(faceMap[face]);
			}
		}

		// --- Texture handling ---
		bool hasNormalMap = !normalMap.empty();
		bool hasHeightMap = !heightMap.empty() && strength > 0.0f;

		Texture* normalTex = new Texture();
		Texture* heightTex = new Texture();

		if (hasNormalMap)
			hasNormalMap = normalTex->load(normalMap);
		if (hasHeightMap)
			hasHeightMap = heightTex->load(heightMap);

		// --- Vertex positions ---
		Vec3 half = size * 0.5f;

		Vec3 p[8] = {
			center + Vec3(-half.x, -half.y, -(-half.z)),
			center + Vec3(half.x, -half.y, -(-half.z)),
			center + Vec3(half.x, half.y, -(-half.z)),
			center + Vec3(-half.x, half.y, -(-half.z)),
			center + Vec3(-half.x, -half.y, -(half.z)),
			center + Vec3(half.x, -half.y, -(half.z)),
			center + Vec3(half.x, half.y, -(half.z)),
			center + Vec3(-half.x, half.y, -(half.z))
		};

		Vec3 normals[6] = {
			Vec3(0, 0, 1),
			Vec3(0, 0, -1),
			Vec3(-1, 0, 0),
			Vec3(1, 0, 0),
			Vec3(0, 1, 0),
			Vec3(0, -1, 0)
		};

		// --- Build geometry ---
		for (int f = 0; f < 6; ++f)
		{
			if (removedFaces.count(f)) continue;

			if (f == 4) // Top face (partitioned)
			{
				int xParts = max(1, (int)partition.x);
				int zParts = max(1, (int)partition.z);

				Vec3 v0 = p[faces[f][0]];
				Vec3 v1 = p[faces[f][1]];
				Vec3 v2 = p[faces[f][2]];
				Vec3 v3 = p[faces[f][3]];

				Frame frame;
				frame.fromVector(normals[f]);

				for (int x = 0; x < xParts; ++x)
				{
					for (int z = 0; z < zParts; ++z)
					{
						float fx0 = (float)x / xParts;
						float fx1 = (float)(x + 1) / xParts;
						float fz0 = (float)z / zParts;
						float fz1 = (float)(z + 1) / zParts;

						Vec3 p00 = lerp(lerp(v0, v1, fx0), lerp(v3, v2, fx0), fz0);
						Vec3 p10 = lerp(lerp(v0, v1, fx1), lerp(v3, v2, fx1), fz0);
						Vec3 p11 = lerp(lerp(v0, v1, fx1), lerp(v3, v2, fx1), fz1);
						Vec3 p01 = lerp(lerp(v0, v1, fx0), lerp(v3, v2, fx0), fz1);

						Vertex verts[4] = {
							{ p00, normals[f], fx0, fz0 },
							{ p10, normals[f], fx1, fz0 },
							{ p11, normals[f], fx1, fz1 },
							{ p01, normals[f], fx0, fz1 }
						};

						if (hasNormalMap)
						{
							for (int i = 0; i < 4; ++i)
							{
								Vec3 sampled = colToVec(normalTex->sample(verts[i].u, verts[i].v));
								Vec3 normal = (sampled * 2.0f - Vec3(1.0f, 1.0f, 1.0f)).normalize();
								std::swap(normal.x, normal.y);
								verts[i].normal = frame.toWorld(normal);
							}
						}

						if (hasHeightMap)
						{
							int gridX[4] = { x, x + 1, x + 1, x };
							int gridZ[4] = { z, z, z + 1, z + 1 };
							for (int i = 0; i < 4; ++i)
							{
								bool isEdge = (gridX[i] == 0 || gridX[i] == xParts) ||
									(gridZ[i] == 0 || gridZ[i] == zParts);
								if (!isEdge)
								{
									float height = heightTex->sample(verts[i].u, verts[i].v).r;
									height -= height * 0.5f;
									verts[i].p = verts[i].p + verts[i].normal * height * strength;
								}
							}
						}

						Triangle t1, t2;
						t1.init(verts[0], verts[1], verts[3], matIndex);
						t2.init(verts[1], verts[2], verts[3], matIndex);

						triangles.push_back(t1);
						triangles.push_back(t2);
					}
				}
			}
			else // Regular face
			{
				Vertex v0 = { p[faces[f][0]], normals[f], 0.0f, 0.0f };
				Vertex v1 = { p[faces[f][1]], normals[f], 1.0f, 0.0f };
				Vertex v2 = { p[faces[f][2]], normals[f], 1.0f, 1.0f };
				Vertex v3 = { p[faces[f][3]], normals[f], 0.0f, 1.0f };

				Triangle t1, t2;
				t1.init(v0, v1, v3, matIndex);
				t2.init(v1, v2, v3, matIndex);

				triangles.push_back(t1);
				triangles.push_back(t2);
			}
		}

		delete normalTex;
		delete heightTex;
	}
};


class PrismMesh : public CustomMesh
{
public:
	PrismMesh(Vec3 center, float radius, float height, unsigned int matIndex = 0)
	{
		float h = height * 0.5f;

		// Flip Z for positions
		Vec3 p0 = center + Vec3(radius, -h, 0);
		Vec3 p1 = center + Vec3(-radius * 0.5f, -h, -radius * sqrtf(3) * 0.5f);
		Vec3 p2 = center + Vec3(-radius * 0.5f, -h, radius * sqrtf(3) * 0.5f);

		Vec3 p3 = p0 + Vec3(0, height, 0);
		Vec3 p4 = p1 + Vec3(0, height, 0);
		Vec3 p5 = p2 + Vec3(0, height, 0);

		// Bottom face (-Y)
		{
			Vec3 normal = Vec3(0, -1, 0);
			Triangle t;
			t.init(Vertex{ p0, normal, 0, 0 }, Vertex{ p2, normal, 1, 0 }, Vertex{ p1, normal, 0.5f, 1 }, matIndex);
			triangles.push_back(t);
		}

		// Top face (+Y) (reversed winding)
		{
			Vec3 normal = Vec3(0, 1, 0);
			Triangle t;
			t.init(Vertex{ p3, normal, 0, 0 }, Vertex{ p4, normal, 1, 0 }, Vertex{ p5, normal, 0.5f, 1 }, matIndex);
			triangles.push_back(t);
		}

		auto addSide = [&](Vec3 a, Vec3 b, Vec3 c, Vec3 d)
			{
				Vec3 normal = (b - a).cross(c - a).normalize();
				Triangle t1, t2;
				// Reverse winding order here too
				t1.init(Vertex{ a, normal, 0, 0 }, Vertex{ c, normal, 1, 1 }, Vertex{ b, normal, 1, 0 }, matIndex);
				t2.init(Vertex{ a, normal, 0, 0 }, Vertex{ d, normal, 0, 1 }, Vertex{ c, normal, 1, 1 }, matIndex);
				triangles.push_back(t1);
				triangles.push_back(t2);
			};

		addSide(p0, p1, p4, p3);
		addSide(p1, p2, p5, p4);
		addSide(p2, p0, p3, p5);
	}
};

class SphereMesh :public CustomMesh
{
public:
	SphereMesh(Vec3 center, float radius, unsigned int segments = 16, unsigned int rings = 16, unsigned int matIndex = 0)
	{
		for (unsigned int y = 0; y < rings; ++y)
		{
			float v0 = (float)y / rings;
			float v1 = (float)(y + 1) / rings;
			float phi0 = v0 * M_PI;
			float phi1 = v1 * M_PI;

			for (unsigned int x = 0; x < segments; ++x)
			{
				float u0 = (float)x / segments;
				float u1 = (float)(x + 1) / segments;

				float theta0 = u0 * 2.0f * M_PI;
				float theta1 = u1 * 2.0f * M_PI;

				Vec3 p00 = center + Vec3(
					radius * sinf(phi0) * cosf(theta0),
					radius * cosf(phi0),
					-radius * sinf(phi0) * sinf(theta0)  // Flip Z here
				);
				Vec3 p01 = center + Vec3(
					radius * sinf(phi1) * cosf(theta0),
					radius * cosf(phi1),
					-radius * sinf(phi1) * sinf(theta0)
				);
				Vec3 p10 = center + Vec3(
					radius * sinf(phi0) * cosf(theta1),
					radius * cosf(phi0),
					-radius * sinf(phi0) * sinf(theta1)
				);
				Vec3 p11 = center + Vec3(
					radius * sinf(phi1) * cosf(theta1),
					radius * cosf(phi1),
					-radius * sinf(phi1) * sinf(theta1)
				);

				Vec3 n00 = (p00 - center).normalize();
				Vec3 n01 = (p01 - center).normalize();
				Vec3 n10 = (p10 - center).normalize();
				Vec3 n11 = (p11 - center).normalize();

				Triangle t1, t2;
				t1.init(Vertex{ p00, n00, u0, v0 }, Vertex{ p01, n01, u0, v1 }, Vertex{ p10, n10, u1, v0 }, matIndex);
				t2.init(Vertex{ p10, n10, u1, v0 }, Vertex{ p01, n01, u0, v1 }, Vertex{ p11, n11, u1, v1 }, matIndex);

				triangles.push_back(t1);
				triangles.push_back(t2);
			}
		}
	}
};
