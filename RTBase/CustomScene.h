#pragma once
#include "Scene.h"
#include "MeshObjects.h"
#include <string.h>
#include "ObjLoader.h"

class CustomScene
{
	unsigned int width;
	unsigned int height;
	float fov;

public:
	const float scale = 3.f; // Scale factor for the scene

	// Constructor with default parameters
	CustomScene(unsigned int width = 1920, unsigned int height = 1080, unsigned int fov = 45.0f)
		: width(width), height(height), fov(fov)
	{
	}

	Scene* createScene(Camera& viewCamera) const
	{
		// create scene
		Scene* scene = new Scene();

		// set up camera
		Matrix P = Matrix::perspective(0.001f, 10000.0f, (float)width / (float)height, fov);

		// set up view matrix
		Vec3 from = Vec3(0, 3.4f, 3.4f) * scale;
		Vec3 to = Vec3(0, 0, 0) * scale;
		Vec3 up(0, 1, 0);
		Matrix V = Matrix::lookAt(from, to, up);
		V = V.invert();

		// initialize camera
		scene->cameraMobile.init(P, width, height);
		scene->cameraMobile.updateView(V);

		// set up view camera
		viewCamera.from = from;
		viewCamera.to = to;
		viewCamera.up = up;
		viewCamera.camera = &scene->cameraMobile;

		// initialize scene data
		std::vector<Triangle> meshTriangles;
		std::vector<BSDF*> meshMaterials;
		std::vector<unsigned int> lights;

		// add materials
		// 0 - plain color
		Texture* tex1 = new Texture(COLOR_LIGHT_GREY);
		meshMaterials.push_back(new PlanktonBSDF(tex1));

		Light* background;
		background = new BackgroundColour(COLOR_BLACK);

		std::vector<Triangle> obj;

		if (loadOBJ("data/Plankton.obj", obj, 0, Vec3(0, 0, 0) * scale, scale))
			meshTriangles.insert(meshTriangles.end(), obj.begin(), obj.end());

		// add sphere
		// SphereMesh sphereMesh(Vec3(0, 0, 0) * scale, 1.0f * scale);
		// meshTriangles.insert(meshTriangles.end(), sphereMesh.triangles.begin(), sphereMesh.triangles.end());

		// add plane mesh
		// PlaneMesh planeMesh(Vec3(0, 0, 0) * scale, Vec3(scale, scale), V3_ONE, 0);
		// meshTriangles.insert(meshTriangles.end(), planeMesh.triangles.begin(), planeMesh.triangles.end());

		// add triangles to scene
		scene->init(meshTriangles, meshMaterials, background);

		// set view camera speed based on scene bounds
		viewCamera.movespeed = (scene->bounds.max - scene->bounds.min).length() * 0.05f;

		// build the scene
		scene->build();

		// set scene bounds
		Vec3 center = (scene->bounds.max + scene->bounds.min) * 0.5f;
		float radius = (scene->bounds.max - use<SceneBounds>().sceneCentre).length();
		use<SceneBounds>().init(center, radius);

		std::cout << "Scene Center: " << center << ", Scene Radius: " << radius << std::endl;

		from = center + V3_BACK * radius * 3;
		to = center;

		V = Matrix::lookAt(from, to, up);
		V = V.invert();

		// initialize camera
		scene->cameraFixed.init(P, width, height);
		scene->cameraFixed.updateView(V);

		return scene;
	}

	bool loadOBJ(std::string path, std::vector<Triangle>& obj, unsigned int material, Vec3 offset = V3_ZERO, float scale = 1.0f) const
	{
		std::vector<OBJVertex> vertices;
		std::vector<OBJTriangle> triangles;

		if (OBJLoader::Load(path, vertices, triangles))
		{
			for (const auto& tri : triangles)
			{
				Vec3 p1 = Vec3(vertices[tri.v1].p.x, vertices[tri.v1].p.y, vertices[tri.v1].p.z);
				Vec3 p2 = Vec3(vertices[tri.v2].p.x, vertices[tri.v2].p.y, vertices[tri.v2].p.z);
				Vec3 p3 = Vec3(vertices[tri.v3].p.x, vertices[tri.v3].p.y, vertices[tri.v3].p.z);

				Vec3 n1 = Vec3(vertices[tri.v1].n.x, vertices[tri.v1].n.y, vertices[tri.v1].n.z);
				Vec3 n2 = Vec3(vertices[tri.v2].n.x, vertices[tri.v2].n.y, vertices[tri.v2].n.z);
				Vec3 n3 = Vec3(vertices[tri.v3].n.x, vertices[tri.v3].n.y, vertices[tri.v3].n.z);

				Vertex v1{ p1 * scale + offset, n1, 0, 0 };
				Vertex v2{ p2 * scale + offset, n2, 1, 0 };
				Vertex v3{ p3 * scale + offset, n3, 0.5f, 1 };

				// Create triangle and add to mesh
				obj.emplace_back(Triangle(v1, v2, v3, material));
			}

			return true;
		}
		else
			std::cerr << "Failed to load OBJ file: " << path << "\n";

		return false;
	}
};