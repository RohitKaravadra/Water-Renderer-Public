#pragma once

enum DRAW_MODE
{
	DM_TRIANGLES,
	DM_NORMALS,
	DM_ALBEDO,
	DM_DIRECT,
	DM_LIGHTS
};

struct SETTINGS
{
	DRAW_MODE drawMode;
	TONEMAP toneMap;
	IMAGE_FILTER filter;

	bool render;					// render enabled
	bool saveRenders;				// save renders to disk

	bool useMultithreading;			// multithreading enabled
	bool adaptiveSampling;			// tile based adaptive sampling
	bool debug;
	bool denoise;					// set to true to denoise the image
	bool printStats;				// print stats to console
	bool debugVMF;					// debug vMF distribution

	unsigned int numThreads;		// number of threads for multithreading
	unsigned int maxBounces;		// max number of bounces for path tracing

	unsigned int initSPP;			// initial samples per pixel
	unsigned int totalSPP;			// total samples per pixel
	float lineWidth;				// line width for triangle edges

	SETTINGS()
	{
		toneMap = TM_NONE;
		filter = FT_BOX;

		render = false;
		saveRenders = false;

		useMultithreading = false;
		adaptiveSampling = false;
		debug = false;
		denoise = false;
		printStats = true;

		numThreads = 3;
		maxBounces = 5;

		initSPP = 10;
		totalSPP = 8192;

		lineWidth = 0.001f;
	}
};

