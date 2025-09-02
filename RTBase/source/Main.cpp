

#include "GEMLoader.h"
#include "Renderer.h"
#include "SceneLoader.h"
#define NOMINMAX
#include "GamesEngineeringBase.h"
#include <unordered_map>
#include "SceneManager.h"

// create settings
SETTINGS createSettings()
{
	SETTINGS settings;

	settings.drawMode = DM_TRIANGLES;
	settings.toneMap = TM_REINHARD_GLOBAL;
	settings.filter = FT_BOX;

	settings.render = false;
	settings.saveRenders = true;
	settings.printStats = false;

	settings.debug = false;
	settings.denoise = true;

	settings.useMultithreading = true;

	settings.adaptiveSampling = false;
	settings.initSPP = 50;
	settings.totalSPP = 500;

	settings.numThreads = 8;
	settings.maxBounces = 6;

	settings.lineWidth = 0.001f;

	return settings;
}

static void saveRender(Renderer& rt, const std::string& type)
{
	if (!rt.settings.saveRenders)
		return;

	const std::string filename = "Image";

	const wchar_t* rendersFolder = L"Renders";
	CreateDirectory(L"Renders", NULL);

	std::string renderFolderStr = std::string(rendersFolder, rendersFolder + wcslen(rendersFolder));
	std::string filepath = renderFolderStr + "/" + filename + type;
	rt.savePNG(filepath);

	std::cout << "Saved render to " << filepath << "\n";
}

class Stats
{
	float totalTime;
	float lastInputTime;
	float renderStartTime;
	float renderTime;
	float deltaTime;

	bool completed;	// flag to check if the render is completed
public:
	bool forceComlete; // flag to force completion for testing

	Stats()
	{
		totalTime = 0;
		lastInputTime = 0;
		renderStartTime = 0;
		renderTime = 0;
		deltaTime = 0;
		completed = false;
		forceComlete = false;
	}

	void reset()
	{
		renderStartTime = totalTime; // reset render start time
		if (completed)
		{
			std::cout << "\n\n\n\n\n";
			completed = false;	// reset completed flag
			forceComlete = false;
		}
	}

	void update(float dt)
	{
		deltaTime = dt;
		totalTime += deltaTime;
	}

	void updateLastInputTime() {
		lastInputTime = totalTime; // reset timer for last input
	}

	void onCompletion()
	{
		if (!completed)
		{
			completed = true;
			renderTime = totalTime - renderStartTime;
		}
	}

	bool canInput() const { return totalTime - lastInputTime > 0.1f; }

	bool isCompleted() const { return completed; }

	void print(Renderer& rt) const
	{
		if (!rt.settings.printStats)
			return;

		float finalTime = completed ? renderTime : totalTime - renderStartTime;
		float progress = rt.settings.render ? rt.getSPP() * 100 / rt.settings.totalSPP : 0;

		// Write stats to console
		std::cout << "\033[F\033[F\033[F\033[F\033[F";
		std::cout << "Progress   : " << progress << "%                                                \n";
		std::cout << "Samples    : " << rt.getSPP() << "                                              \n";
		std::cout << "Time       : " << deltaTime << "                                                \n";
		std::cout << "FPS        : " << (deltaTime > 0 ? 1.0f / deltaTime : FLT_MAX) << "             \n";
		std::cout << "Total time : " << std::roundf(totalTime) << " sec                               \n";
	}
};

int main()
{
	// Load scene
	SceneManager sceneManager;
	sceneManager.load(SCENES::CUSTOM);

	// Create canvas
	GamesEngineeringBase::Window canvas;
	canvas.create((unsigned int)sceneManager.curScene->camera.width, (unsigned int)sceneManager.curScene->camera.height, "Tracer", 1.0f);

	// Create ray tracer
	Renderer rt;
	rt.init(sceneManager.curScene, &canvas, createSettings());

	// Create timer
	GamesEngineeringBase::Timer timer;
	Stats stats;

	bool running = true;
	AOV aov;

	if (rt.settings.printStats)
		std::cout << "\n\n\n\n\n";

	while (running)
	{
		stats.update(timer.dt());

		canvas.checkInput(); // Check for input

		// Check if the user wants to quit
		if (canvas.isQuitRequested() || canvas.keyPressed(VK_ESCAPE))
		{
			running = false;
			continue;
		}

		// Update camera and check if it has changed (reset if it has)
		if (sceneManager.viewcamera->update(canvas))
		{
			rt.clear();
			stats.reset();
		}

		if (stats.canInput())
		{
			bool inputChanged = true;

			if (canvas.keyPressed('R'))
				rt.toggleRender();
			else if (canvas.keyPressed('V'))
				rt.toggleLightVMF();
			else if (canvas.keyPressed(VK_SPACE))
				rt.cycleDrawMode();
			else if (canvas.keyPressed(VK_BACK))
				stats.forceComlete = true;
			else
				inputChanged = false;

			if (inputChanged)
			{
				stats.reset();
				stats.updateLastInputTime();
			}
		}

		canvas.clear();

		if (!stats.isCompleted())
		{
			if (rt.settings.useMultithreading)
				rt.renderMT();
			else
				rt.render();

			stats.print(rt);

			if (stats.forceComlete || rt.settings.render && rt.settings.totalSPP <= rt.getSPP())
			{
				stats.onCompletion();

				rt.draw();
				saveRender(rt, "_render.png");

				// denoising
				if (rt.settings.denoise)
				{
					rt.createAOV(aov);
					Denoiser denoiser(aov.width, aov.height);
					denoiser.denoise(aov);
				}

				rt.draw(aov);
				saveRender(rt, "_denoised.png");
			}
			rt.draw();
		}

		else
		{
			if (rt.settings.denoise)
				rt.draw(aov);
			else
				rt.draw();
		}

		// debug for camera position
		if (!rt.settings.printStats && canvas.keyPressed('P'))
		{
			std::cout << "Camera Position : " << sceneManager.viewcamera->from << "\n";
			std::cout << "Camera Target   : " << sceneManager.viewcamera->to << "\n\n";
		}

		// draw the image
		canvas.present();
	}

	return 0;
}