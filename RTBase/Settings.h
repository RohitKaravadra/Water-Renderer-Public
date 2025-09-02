#pragma once

#include "Imaging.h"

enum DRAW_MODE
{
	DM_VMF,
	DM_TRIANGLES,
	DM_NORMALS,
	DM_ALBEDO,
};

struct SETTINGS
{
	DRAW_MODE drawMode;
	TONEMAP toneMap;
	IMAGE_FILTER filter;

	bool render;					// render enabled
	bool paused;					// paused
	bool useMultithreading;			// multithreading enabled
	bool printStats;				// print stats to console

	unsigned int numThreads;		// number of threads for multithreading
	unsigned int maxTraces;		// max number of bounces for path tracing
	unsigned int tracePerFrame;		// number of ray traces per frame	

	unsigned int totalSPP;			// total samples per pixel
	float lineWidth;				// line width for triangle edges

	SETTINGS()
	{
		drawMode = DM_TRIANGLES;
		toneMap = TM_NONE;
		filter = FT_BOX;

		render = false;
		paused = false;
		useMultithreading = false;
		printStats = true;
		tracePerFrame = 1;

		numThreads = 3;
		maxTraces = 5;

		totalSPP = 8192;

		lineWidth = 0.001f;
	}
};

