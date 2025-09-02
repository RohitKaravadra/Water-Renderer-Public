#pragma once

#include <mutex>
#include <algorithm>
#include <cmath>

#include "Core.h"
#include "Imaging.h"
#include "Scene.h"

struct DebugLine
{
	Vec3 start, end;
	Color color1, color2;
	bool drawEnd, blend;

	DebugLine() : start(0, 0, 0), end(0, 0, 0), color1(COLOR_WHITE), color2(COLOR_WHITE), drawEnd(false), blend(false) {}
	DebugLine(Vec3 s, Vec3 e, Color c1 = COLOR_CYAN, bool drawEnd = false, Color c2 = COLOR_BLACK, bool blend = false)
		: start(s), end(e), color1(c1), color2(c2), drawEnd(drawEnd), blend(blend) {
	}
};

class Debugger
{
	Film* film;                             // reference to the film for rendering
	Scene* scene;                           // reference to the scene for camera info

	std::vector<DebugLine> debugLines;      // stores the debug lines for rendering
	std::mutex debugMutex;                  // mutex for thread safety

	void clipToScreen(Vec3& p1, Vec3& p2, int width, int height)
	{
		// Define region codes
		const int INSIDE = 0;
		const int LEFT = 1;
		const int RIGHT = 2;
		const int BOTTOM = 4;
		const int TOP = 8;
		const int _NEAR = 16; // z < 0 (behind camera)

		// Function to compute region code
		auto computeCode = [&](const Vec3& p) -> int {
			int code = INSIDE;

			if (p.x < 0) code |= LEFT;
			else if (p.x > width) code |= RIGHT;

			if (p.y < 0) code |= BOTTOM;
			else if (p.y > height) code |= TOP;

			if (p.z < 0) code |= _NEAR; // Behind camera

			return code;
			};

		int code1 = computeCode(p1);
		int code2 = computeCode(p2);

		int iteration = 0;
		const int MAX_ITERATIONS = 10;

		while (iteration++ < MAX_ITERATIONS) {
			if ((code1 | code2) == 0) {
				// Both points inside
				break;
			}
			if (code1 & code2) {
				// Both points outside the same edge
				break;
			}

			int codeOut = code1 ? code1 : code2;
			Vec3* pOut = code1 ? &p1 : &p2;

			float t = 0.0f;
			Vec3 intersection;

			// Calculate intersection with the appropriate plane
			if (codeOut & TOP) {
				t = (height - p1.y) / (p2.y - p1.y);
				intersection.x = p1.x + (p2.x - p1.x) * t;
				intersection.y = height;
				intersection.z = p1.z + (p2.z - p1.z) * t;
			}
			else if (codeOut & BOTTOM) {
				t = (0 - p1.y) / (p2.y - p1.y);
				intersection.x = p1.x + (p2.x - p1.x) * t;
				intersection.y = 0;
				intersection.z = p1.z + (p2.z - p1.z) * t;
			}
			else if (codeOut & RIGHT) {
				t = (width - p1.x) / (p2.x - p1.x);
				intersection.x = width;
				intersection.y = p1.y + (p2.y - p1.y) * t;
				intersection.z = p1.z + (p2.z - p1.z) * t;
			}
			else if (codeOut & LEFT) {
				t = (0 - p1.x) / (p2.x - p1.x);
				intersection.x = 0;
				intersection.y = p1.y + (p2.y - p1.y) * t;
				intersection.z = p1.z + (p2.z - p1.z) * t;
			}
			else if (codeOut & _NEAR) {
				// Intersection with z=0 plane (camera plane)
				t = (0 - p1.z) / (p2.z - p1.z);
				intersection.x = p1.x + (p2.x - p1.x) * t;
				intersection.y = p1.y + (p2.y - p1.y) * t;
				intersection.z = 0;
			}

			// Replace the outside point with the intersection point
			*pOut = intersection;

			// Update the region code for the modified point
			if (codeOut == code1) {
				code1 = computeCode(p1);
			}
			else {
				code2 = computeCode(p2);
			}
		}
	}

	void drawLine(const DebugLine& line)
	{
		Vec3 p0 = line.start;
		Vec3 p1 = line.end;

		// Project both points (modifies them in-place)
		bool startInside = scene->camera->projectOntoCamera(p0);
		bool endInside = scene->camera->projectOntoCamera(p1);

		// If both points are behind the camera (z < 0), skip entirely
		if (p0.z < 0 || p1.z < 0) {
			return;
		}

		// Convert to screen coordinates
		// p0 and p1 are now in screen coordinates with z = distance from camera

		// Clip to screen boundaries AND near plane (z=0)
		clipToScreen(p0, p1, film->width, film->height);

		// Check if the line is completely outside after clipping
		if ((p0.x < 0 && p1.x < 0) || (p0.x > film->width && p1.x > film->width) ||
			(p0.y < 0 && p1.y < 0) || (p0.y > film->height && p1.y > film->height) ||
			(p0.z < 0 && p1.z < 0)) {
			return;
		}

		// If either point is behind camera after clipping, skip
		if (p0.z < 0 || p1.z < 0) {
			return;
		}

		// Clamp points to screen boundaries
		p0.x = clamp(p0.x, 0.0f, static_cast<float>(film->width - 1));
		p0.y = clamp(p0.y, 0.0f, static_cast<float>(film->height - 1));
		p1.x = clamp(p1.x, 0.0f, static_cast<float>(film->width - 1));
		p1.y = clamp(p1.y, 0.0f, static_cast<float>(film->height - 1));

		p0.z = p1.z = 0; // We don't need depth for 2D line drawing
		float totalDist = (p1 - p0).length();

		int x0 = static_cast<int>(std::round(p0.x));
		int y0 = static_cast<int>(std::round(p0.y));
		int x1 = static_cast<int>(std::round(p1.x));
		int y1 = static_cast<int>(std::round(p1.y));

		if (x0 == x1 && y0 == y1) {
			return;
		}

		// Bresenham line drawing
		int dx = std::abs(x1 - x0);
		int dy = std::abs(y1 - y0);
		int sx = (x0 < x1) ? 1 : -1;
		int sy = (y0 < y1) ? 1 : -1;
		int err = dx - dy;

		int x = x0;
		int y = y0;

		while (true) {
			if (x >= 0 && x < film->width && y >= 0 && y < film->height)
			{
				Color col = line.color1;
				if (line.blend)
				{
					float currentDist = (Vec3(x, y) - p0).length();
					float t = clamp(currentDist / totalDist, 0.0f, 1.0f);
					col = line.color1.blend(line.color2, t);
				}
				film->splat(x, y, col);
			}

			if (x == x1 && y == y1) break;

			int e2 = 2 * err;
			if (e2 > -dy) {
				err -= dy;
				x += sx;
			}
			if (e2 < dx) {
				err += dx;
				y += sy;
			}
		}
	}

	void drawPoint(Vec3 point, Color color = COLOR_RED)
	{
		if (scene->camera->projectOntoCamera(point))
		{
			if (point.z < 0)
				return; // Point is behind the camera

			int x = static_cast<int>(std::round(point.x));
			int y = static_cast<int>(std::round(point.y));

			// Draw a small cross for the point
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					if (std::abs(i) + std::abs(j) <= 1) { // Cross pattern
						int px = x + i;
						int py = y + j;
						if (px >= 0 && px < film->width && py >= 0 && py < film->height) {
							film->splat(px, py, color);
						}
					}
				}
			}
		}
	}

public:
	void init(Film*& _film, Scene*& _scene)
	{
		film = _film;
		scene = _scene;
	}

	void addLine(const DebugLine& debugLine)
	{
		std::lock_guard<std::mutex> lock(debugMutex);
		debugLines.push_back(debugLine);
	}

	void clear() {
		debugLines.clear();
	}

	void draw()
	{
		for (const DebugLine& line : debugLines)
		{
			drawLine(line);
			if (line.drawEnd)
			{
				drawPoint(line.end, COLOR_YELLOW);
				drawPoint(line.start, COLOR_YELLOW);
			}
		}
	}
};