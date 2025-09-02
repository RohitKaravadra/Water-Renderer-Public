#pragma once

#include <string.h>

#include "Scene.h"
#include "MeshObjects.h"
#include "ObjLoader.h"

class CustomScene
{
	unsigned int width;
	unsigned int height;
	float fov;

	enum SceneTypes
	{
		SCENE_POOL,
		SCENE_POND,
		SCENE_OCEAN,
		SCENE_SPHERE,
		SCENE_WATERCUBE
	};

	enum CameraViews
	{
		CSC_DEFAULT,
		CSC_UNDERWATER
	};

public:
	const bool usePond = true;
	const bool useEnvMap = true;
	const bool useSphere = false;

	SceneTypes sceneType = SCENE_OCEAN;
	// CameraViews cameraView = CSC_DEFAULT;
	CameraViews cameraView = CSC_UNDERWATER;

	const float scale = 1.f; // Scale factor for the scene

	// Constructor with default parameters
	CustomScene(unsigned int width = 1920, unsigned int height = 1080, unsigned int fov = 45.0f)
		: width(width), height(height), fov(fov)
	{
	}

	void createPondScene(Scene*& scene, Camera& viewCamera, const Matrix& projMat) const
	{
		// set up view matrix
		Vec3 from, to;
		if (cameraView == CSC_UNDERWATER)
		{
			from = Vec3(-4.28039f, 0.3933f, 0.737121f) * scale;
			to = Vec3(0, 0, 0) * scale;
		}
		else
		{
			from = Vec3(-5.31867, 4.14358, 0.931809) * scale;
			to = Vec3(-4.62394, 3.43433, 0.81217) * scale;
		}

		Vec3 up(0, 1, 0);
		Matrix V = Matrix::lookAt(from, to, up);
		V = V.invert();

		// initialize camera
		scene->camera.init(projMat, width, height);
		scene->camera.updateView(V);

		// set up view camera
		viewCamera.from = from;
		viewCamera.to = to;
		viewCamera.up = up;
		viewCamera.camera = &scene->camera;

		// initialize scene data
		std::vector<Triangle> sceneTriangles;
		std::vector<BSDF*> sceneMaterials;
		std::vector<unsigned int> lights;

		// add materials
		// 0 - plain color
		Texture* tex1 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new DiffuseBSDF(tex1));

		// 1 - water
		Texture* tex2 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new WaterBSDF(tex2));

		Texture* env = new Texture();
		env->load("data/Sky_3.hdr");
		Light* background = new EnvironmentMap(env);

		CubeMesh volumeMesh(Vec3(0, 0.25f, 0) * scale, Vec3(15, 0.6f, 5) * scale, V3_ONE * 200,
			"data/water_n.png", "data/water_h.png", 0.1f * scale, 1);
		scene->addVolume(volumeMesh.triangles);

		// load pond model
		std::vector<Triangle> meshTriangles;
		if (loadOBJ("data/Pond.obj", meshTriangles, 0))
			sceneTriangles.insert(sceneTriangles.end(), meshTriangles.begin(), meshTriangles.end());

		// add materials to scene
		scene->init(sceneTriangles, sceneMaterials, background);
	}

	void createPoolScene(Scene*& scene, Camera& viewCamera, const Matrix& projMat) const
	{
		Vec3 from, to;
		if (cameraView == CSC_UNDERWATER)
		{
			from = Vec3(0.736191, -0.401973, 3.09712) * scale;
			to = Vec3(0.823338, -0.414424, 2.101) * scale;
		}
		else
		{
			from = Vec3(-1.20738f, 3.49964f, 5.53129f) * scale;
			to = Vec3(-0.913338f, 2.9889f, 4.72341f) * scale;
		}

		Vec3 up(0, 1, 0);
		Matrix V = Matrix::lookAt(from, to, up);
		V = V.invert();

		// initialize camera
		scene->camera.init(projMat, width, height);
		scene->camera.updateView(V);

		// set up view camera
		viewCamera.from = from;
		viewCamera.to = to;
		viewCamera.up = up;
		viewCamera.camera = &scene->camera;

		// initialize scene data
		std::vector<Triangle> sceneTriangles;
		std::vector<BSDF*> sceneMaterials;
		std::vector<unsigned int> lights;

		// add materials
		// 0 - plain color
		Texture* tex1 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new DiffuseBSDF(tex1));

		// 1 - water
		Texture* tex2 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new WaterBSDF(tex2));

		Texture* env = new Texture();
		env->load("data/Sky_3.hdr");
		Light* background = new EnvironmentMap(env);


		CubeMesh volumeMesh(Vec3(0.5, -1.f, 0.f) * scale, Vec3(5, 3, 10) * scale, V3_ONE * 50,
			"data/water_n.png", "data/water_h.png", 0.1f * scale, 1);
		scene->addVolume(volumeMesh.triangles);

		// load pond model
		std::vector<Triangle> meshTriangles;
		if (loadOBJ("data/Pool.obj", meshTriangles, 0, Vec3(0, 0), 0.2f))
			sceneTriangles.insert(sceneTriangles.end(), meshTriangles.begin(), meshTriangles.end());

		// add materials to scene
		scene->init(sceneTriangles, sceneMaterials, background);
	}

	void createOceanScene(Scene*& scene, Camera& viewCamera, const Matrix& projMat) const
	{
		// set up view matrix

		Vec3 from, to;
		if (cameraView == CSC_UNDERWATER)
		{
			from = Vec3(-22.697, -0.74157, 19.0235) * scale;
			to = Vec3(-22.9478, -0.494846, 19.9596) * scale;
		}
		else
		{
			from = Vec3(-25.5116, 22.6373, 30.058) * scale;
			to = Vec3(-25.0244, 21.985, 29.4774) * scale;
		}


		Vec3 up(0, 1, 0);
		Matrix V = Matrix::lookAt(from, to, up);
		V = V.invert();

		// initialize camera
		scene->camera.init(projMat, width, height);
		scene->camera.updateView(V);

		// set up view camera
		viewCamera.from = from;
		viewCamera.to = to;
		viewCamera.up = up;
		viewCamera.camera = &scene->camera;

		// initialize scene data
		std::vector<Triangle> sceneTriangles;
		std::vector<BSDF*> sceneMaterials;
		std::vector<unsigned int> lights;

		// 0 - plain color
		Texture* tex1 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new DiffuseBSDF(tex1));

		// 1 - water
		Texture* tex2 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new WaterBSDF(tex2));

		Texture* env = new Texture();
		env->load("data/Sky_3.hdr");
		Light* background = new EnvironmentMap(env);

		CubeMesh volumeMesh(Vec3(0, 0.f, 0) * scale, Vec3(100, 10, 100) * scale, V3_ONE * 200,
			"data/water_n.png", "data/water_h.png", 0.1f * scale, 1);
		scene->addVolume(volumeMesh.triangles);

		// add a sphere inside the volume near to surface
		SphereMesh sphereMesh(Vec3(-25.f, 2.f, 29.f) * scale, 1.f * scale, 32, 32, 0);
		sceneTriangles.insert(sceneTriangles.end(), sphereMesh.triangles.begin(), sphereMesh.triangles.end());

		// add materials to scene
		scene->init(sceneTriangles, sceneMaterials, background);
	}

	void createSphereScene(Scene*& scene, Camera& viewCamera, const Matrix& projMat) const
	{
		// default
		Vec3 from = Vec3(0, 3.4f, 7.4f) * scale;
		Vec3 to = Vec3(0, 2, 0) * scale;

		Vec3 up(0, 1, 0);
		Matrix V = Matrix::lookAt(from, to, up);
		V = V.invert();

		// initialize camera
		scene->camera.init(projMat, width, height);
		scene->camera.updateView(V);

		// set up view camera
		viewCamera.from = from;
		viewCamera.to = to;
		viewCamera.up = up;
		viewCamera.camera = &scene->camera;

		// initialize scene data
		std::vector<Triangle> sceneTriangles;
		std::vector<BSDF*> sceneMaterials;
		std::vector<unsigned int> lights;

		// add materials
		// 0 - plain color
		Texture* tex1 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new DiffuseBSDF(tex1));

		// 1 - water
		Texture* tex2 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new WaterBSDF(tex2));


		Light* background = new BackgroundColour(COLOR_BLACK);

		// 2 - light material
		BSDF* lightMaterial = new DiffuseBSDF(new Texture(COLOR_WHITE));
		lightMaterial->addLight(COLOR_WHITE * 17.f);
		sceneMaterials.push_back(lightMaterial);

		// add light mesh
		PlaneMesh lightMesh(Vec3(0, 4.99f, 0) * scale, Vec3(2, 2) * scale, V3_ONE, 2);
		lightMesh.invertNormals(); // invert normals for light mesh
		sceneTriangles.insert(sceneTriangles.end(), lightMesh.triangles.begin(), lightMesh.triangles.end());

		SphereMesh sphereMesh(Vec3(0, 2.f, 0) * scale, scale, 32, 32, 1);
		sceneTriangles.insert(sceneTriangles.end(), sphereMesh.triangles.begin(), sphereMesh.triangles.end());

		// add wall on back (cube mesh)
		CubeMesh boundaryBox(Vec3(0, 5.f, 0) * scale, Vec3(15, 10, 10) * scale, V3_ONE, 0, "back");
		boundaryBox.invertNormals(); // invert normals for back wall
		sceneTriangles.insert(sceneTriangles.end(), boundaryBox.triangles.begin(), boundaryBox.triangles.end());

		// add materials to scene
		scene->init(sceneTriangles, sceneMaterials, background);
	}

	void createWaterCubeScene(Scene*& scene, Camera& viewCamera, const Matrix& projMat) const
	{
		// set up view matrix
		Vec3 from = Vec3(0, 3.4f, 7.4f) * scale;
		Vec3 to = Vec3(0, 2, 0) * scale;

		Vec3 up(0, 1, 0);
		Matrix V = Matrix::lookAt(from, to, up);
		V = V.invert();

		// initialize camera
		scene->camera.init(projMat, width, height);
		scene->camera.updateView(V);

		// set up view camera
		viewCamera.from = from;
		viewCamera.to = to;
		viewCamera.up = up;
		viewCamera.camera = &scene->camera;

		// initialize scene data
		std::vector<Triangle> sceneTriangles;
		std::vector<BSDF*> sceneMaterials;
		std::vector<unsigned int> lights;

		// add materials
		// 0 - plain color
		Texture* tex1 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new DiffuseBSDF(tex1));

		// 1 - water
		Texture* tex2 = new Texture(COLOR_WHITE);
		sceneMaterials.push_back(new WaterBSDF(tex2));

		Light* background = new BackgroundColour(COLOR_BLACK);

		// 2 - light material
		BSDF* lightMaterial = new DiffuseBSDF(new Texture(COLOR_WHITE));
		lightMaterial->addLight(COLOR_WHITE * 17.f);
		sceneMaterials.push_back(lightMaterial);

		// add light mesh
		PlaneMesh lightMesh(Vec3(0, 4.99f, 0) * scale, Vec3(2, 2) * scale, V3_ONE, 2);
		lightMesh.invertNormals(); // invert normals for light mesh
		sceneTriangles.insert(sceneTriangles.end(), lightMesh.triangles.begin(), lightMesh.triangles.end());

		// add wall on back (cube mesh)
		CubeMesh boundaryBox(Vec3(0, 5.f, 0) * scale, Vec3(15, 10, 10) * scale, V3_ONE, 0, "back");
		boundaryBox.invertNormals(); // invert normals for back wall
		sceneTriangles.insert(sceneTriangles.end(), boundaryBox.triangles.begin(), boundaryBox.triangles.end());


		CubeMesh volumeMesh(Vec3(0, 2.f, 0) * scale, Vec3(6, 2, 3) * scale, V3_ONE * 20,
			"data/water_n.png", "data/water_h.png", 0.1f * scale, 1);
		sceneTriangles.insert(sceneTriangles.end(), volumeMesh.triangles.begin(), volumeMesh.triangles.end());

		// add materials to scene
		scene->init(sceneTriangles, sceneMaterials, background);
	}

	Scene* createScene(Camera& viewCamera) const
	{
		// create scene
		Scene* scene = new Scene();

		// set up camera
		Matrix projMat = Matrix::perspective(0.001f, 10000.0f, (float)width / (float)height, fov);

		switch (sceneType)
		{
		case SCENE_POOL:
			createPoolScene(scene, viewCamera, projMat);
			break;
		case SCENE_POND:
			createPondScene(scene, viewCamera, projMat);
			break;
		case SCENE_OCEAN:
			createOceanScene(scene, viewCamera, projMat);
			break;
		case SCENE_SPHERE:
			createSphereScene(scene, viewCamera, projMat);
			break;
		case SCENE_WATERCUBE:
			createWaterCubeScene(scene, viewCamera, projMat);
			break;
		default:
			createSphereScene(scene, viewCamera, projMat);
		}

		// set view camera speed based on scene bounds
		viewCamera.movespeed = (scene->bounds.max - scene->bounds.min).length() * 0.05f;

		// build the scene
		scene->build();

		// set scene bounds
		use<SceneBounds>().sceneCentre = (scene->bounds.max + scene->bounds.min) * 0.5f;
		use<SceneBounds>().sceneRadius = (scene->bounds.max - use<SceneBounds>().sceneCentre).length();

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