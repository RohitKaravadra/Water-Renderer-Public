#pragma once
#include <string>
#include "SceneLoader.h"
#include "CustomScene.h"

static enum SCENES
{
	CUSTOM,
	CORNELL_BOX
};

class SceneManager
{
public:

	Scene* curScene = nullptr;
	Camera* viewcamera;

	SceneManager()
	{
		viewcamera = new Camera();
	}

	void unload()
	{
		if (curScene)
		{
			delete curScene;
			curScene = nullptr;
		}
	}

	void load(SCENES scene)
	{
		unload();
		switch (scene)
		{
		case SCENES::CUSTOM:
		{
			std::cout << "Loading scene Custom Scene" << std::endl;
			CustomScene customScene(1080, 720);
			curScene = customScene.createScene(*viewcamera);
		}
		break;
		case SCENES::CORNELL_BOX:
			curScene = loadScene("scenes/cornell-box", *viewcamera);
			break;
		default:
			std::cerr << "Scene not found!" << std::endl;
		}
	}

	~SceneManager()
	{
		unload();
		if (viewcamera)
		{
			delete viewcamera;
			viewcamera = nullptr;
		}
	}
};