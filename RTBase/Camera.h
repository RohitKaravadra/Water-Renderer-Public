#pragma once

#include "Core.h"
#include "Geometry.h"
#include "GamesEngineeringBase.h"

using GamesEngineeringBase::Window;

class SceneCamera
{
public:
	Matrix projectionMatrix;
	Matrix inverseProjectionMatrix;
	Matrix camera;
	Matrix cameraToView;
	float width = 0;
	float height = 0;
	Vec3 origin;
	Vec3 viewDirection;
	float Afilm;
	void init(Matrix ProjectionMatrix, int screenwidth, int screenheight)
	{
		projectionMatrix = ProjectionMatrix;
		inverseProjectionMatrix = ProjectionMatrix.invert();
		width = (float)screenwidth;
		height = (float)screenheight;
		float Wlens = (2.0f / ProjectionMatrix.a[1][1]);
		float aspect = ProjectionMatrix.a[0][0] / ProjectionMatrix.a[1][1];
		float Hlens = Wlens * aspect;
		Afilm = Wlens * Hlens;
	}
	void updateView(Matrix V)
	{
		camera = V;
		cameraToView = V.invert();
		origin = camera.mulPoint(Vec3(0, 0, 0));
		viewDirection = inverseProjectionMatrix.mulPointAndPerspectiveDivide(Vec3(0, 0, 1));
		viewDirection = camera.mulVec(viewDirection);
		viewDirection = viewDirection.normalize();
	}

	/// <summary>
	/// Geberates a ray from camera origin
	/// </summary>
	/// <param name="x"> x coordinate of screen pixel </param>
	/// <param name="y"> y coordinate of screen pixel</param>
	/// <returns> A ray from camera origin with direction to the point on near plane </returns>
	Ray generateRay(float x, float y)
	{
		float xprime = x / width;
		float yprime = 1.0f - (y / height);
		xprime = (xprime * 2.0f) - 1.0f;
		yprime = (yprime * 2.0f) - 1.0f;
		Vec3 dir(xprime, yprime, 1.0f);
		dir = inverseProjectionMatrix.mulPoint(dir);
		dir = camera.mulVec(dir);
		return Ray(origin, dir.normalize());
	}

	bool projectOntoCamera(const Vec3& p, float& x, float& y)
	{
		Vec3 pview = cameraToView.mulPoint(p);
		Vec3 pproj = projectionMatrix.mulPointAndPerspectiveDivide(pview);
		x = (pproj.x + 1.0f) * 0.5f;
		y = (pproj.y + 1.0f) * 0.5f;

		bool inside = !(x < 0 || x > 1.0f || y < 0 || y > 1.0f);

		x = x * width;
		y = 1.0f - y;
		y = y * height;

		return inside;
	}

	bool projectOntoCamera(Vec3& p)
	{
		Vec3 pview = cameraToView.mulPoint(p);
		p = projectionMatrix.mulPointAndPerspectiveDivide(pview);
		p.x = (p.x + 1.0f) * 0.5f;
		p.y = (p.y + 1.0f) * 0.5f;

		bool inside = !(p.x < 0 || p.x > 1.0f || p.y < 0 || p.y > 1.0f);

		p.x = p.x * width;
		p.y = 1.0f - p.y;
		p.y = p.y * height;

		p.z = -pview.z;
		return inside;
	}
};

class Camera
{
public:
	Vec3 from;
	Vec3 to;
	Vec3 up;
	SceneCamera* camera = NULL;
	float movespeed = 1.0f;
	float speedUpdateThreshold = 1.f;
	float rotspeed = 5.0f;
	Camera()
	{
		rotspeed = 5.0f;
	}

	void forward()
	{
		Vec3 dir = (to - from).normalize();
		from = from + dir * movespeed;
		to = from + dir;
	}

	void back()
	{
		Vec3 dir = (to - from).normalize();
		from = from - dir * movespeed;
		to = from + dir;
	}

	void right()
	{
		Vec3 dir = (to - from).normalize();
		Vec3 r = Cross(dir, up);
		from = from + r * movespeed;
		to = from + dir;
	}

	void left()
	{
		Vec3 dir = (to - from).normalize();
		Vec3 r = Cross(dir, up);
		from = from - r * movespeed;
		to = from + dir;
	}

	void flyUp()
	{
		Vec3 dir = up * movespeed;
		from = from + dir;
		to = to + dir;
	}

	void flyDown()
	{
		Vec3 dir = up * movespeed;
		from = from - dir;
		to = to - dir;
	}

	void yaw(float angle)
	{
		Vec3 dir = to - from;
		dir = dir.normalize();

		float rad = angle * (M_PI / 180.0f);
		float cosTheta = cosf(rad);
		float sinTheta = sinf(rad);

		Vec3 k = up.normalize();

		Vec3 rotated = (dir * cosTheta) + (k.cross(dir) * sinTheta) + (k * (k.dot(dir) * (1 - cosTheta)));
		dir = rotated;

		to = from + dir;
	}

	void pitch(float angle)
	{
		Vec3 dir = to - from;
		dir = dir.normalize();

		float rad = angle * (M_PI / 180.0f);
		float cosTheta = cosf(rad);
		float sinTheta = sinf(rad);

		Vec3 k = dir.cross(up).normalize();

		Vec3 rotated = (dir * cosTheta) + (k.cross(dir) * sinTheta) + (k * (k.dot(dir) * (1 - cosTheta)));
		dir = rotated;

		to = from + dir;
	}

	void roll(float angle)
	{
		Vec3 dir = (to - from).normalize();

		float rad = angle * (M_PI / 180.0f);
		float cosTheta = cosf(rad);
		float sinTheta = sinf(rad);

		Vec3 k = dir;

		Vec3 rotated = (up * cosTheta) + (k.cross(up) * sinTheta) + (k * (k.dot(up) * (1 - cosTheta)));
		up = rotated.normalize();
	}
	void updateCamera()
	{
		Matrix V = Matrix::lookAt(from, to, up);
		V = V.invert();
		camera->updateView(V);
	}

	bool update(Window& canvas)
	{
		bool change = false;

		// movement
		if (canvas.keyPressed('W'))
		{
			forward();
			change = true;
		}
		if (canvas.keyPressed('S'))
		{
			back();
			change = true;
		}
		if (canvas.keyPressed('A'))
		{
			left();
			change = true;
		}
		if (canvas.keyPressed('D'))
		{
			right();
			change = true;
		}
		if (canvas.keyPressed('E'))
		{
			flyUp();
			change = true;
		}
		if (canvas.keyPressed('Q'))
		{
			flyDown();
			change = true;
		}

		// rotation
		if (canvas.keyPressed(VK_LEFT))
		{
			if (canvas.keyPressed(VK_SHIFT))
				roll(-rotspeed);
			else
				yaw(rotspeed);
			change = true;
		}
		if (canvas.keyPressed(VK_RIGHT))
		{
			if (canvas.keyPressed(VK_SHIFT))
				roll(rotspeed);
			else
				yaw(-rotspeed);
			change = true;
		}
		if (canvas.keyPressed(VK_UP))
		{
			pitch(rotspeed);
			change = true;
		}
		if (canvas.keyPressed(VK_DOWN))
		{
			pitch(-rotspeed);
			change = true;
		}

		if (change)
			updateCamera();

		int mouseWheel = canvas.getMouseWheel();
		if (mouseWheel != 0)
		{
			movespeed += (mouseWheel < 1 ? -1 : 1) * speedUpdateThreshold;
			if (movespeed < 0.01f)
				movespeed = 0.01f; // Prevent negative or zero speed
		}

		return change;
	}
};