

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

	settings.drawMode = DM_ALBEDO;
	settings.toneMap = TM_REINHARD_GLOBAL;
	settings.filter = FT_BOX;

	settings.printStats = false;

	settings.useMultithreading = true;
	settings.totalSPP = 1000;

	settings.tracePerFrame = 100;
	settings.maxTraces = 20000;

	settings.numThreads = 8;

	settings.lineWidth = 0.001f;

	return settings;
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

	Stats()
	{
		totalTime = 0;
		lastInputTime = 0;
		renderStartTime = 0;
		renderTime = 0;
		deltaTime = 0;
		completed = false;
	}

	void reset()
	{
		renderStartTime = totalTime; // reset render start time
		if (completed)
		{
			std::cout << "\n\n\n\n\n";
			completed = false;	// reset completed flag
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
	canvas.create((unsigned int)sceneManager.curScene->camera->width,
		(unsigned int)sceneManager.curScene->camera->height, "Tracer", 1.0f);

	// Create ray tracer
	Renderer rt;
	rt.init(sceneManager.curScene, &canvas, createSettings());

	// Create timer
	GamesEngineeringBase::Timer timer;
	Stats stats;

	bool running = true;

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

			if (canvas.keyPressed('R'))rt.toggleRender();
			else if (canvas.keyPressed(VK_SPACE))
				rt.cycleDrawMode();
			else if (canvas.keyPressed('P'))
				rt.saveTraceData();
			else if (canvas.keyPressed(VK_TAB))
				rt.settings.paused = !rt.settings.paused;
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

			if (rt.settings.render && rt.settings.totalSPP <= rt.getSPP())
				stats.onCompletion();
		}

		// draw the image
		rt.draw();
		canvas.present();
	}

	return 0;
}