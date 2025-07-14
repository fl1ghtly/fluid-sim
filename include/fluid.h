#pragma once
#include <ctime>
#include <vector>
#include <algorithm>
#include "Vector2f.h"

class Fluid {
	private:
		int width;
		int height;
		int numParticles;
		float radius;
		std::vector<Vector2f> position;
		std::vector<Vector2f> velocity;
		std::vector<float> mass;
		std::vector<float> density;
		std::vector<float> pressure;
		clock_t oldTime;
		void initializeParticleValues();
		void calculateDensity();
		void calculatePressure();
		std::vector<float> calculateAlphaFactors();
		void applyNonPressureForce(float dt);
		void correctDensityError(std::vector<float> alpha, float dt);
		void correctDivergenceError(std::vector<float> alpha, float dt);
		float smoothingKernel(float dist, float smoothingLength);
		Vector2f smoothingGradient(Vector2f r, float smoothingLength);
		template <typename T>
		Vector2f gradient(int particleIndex, std::vector<T> field);
		float divergence(int particleIndex, std::vector<Vector2f> field);
	public:
		Fluid(int width, int heigh, int numParticles, float radius);
		void update();
		std::vector<Vector2f> getPosition();
		void initializeParticleGrid(int x, int y, int spacing, int width, int height);
		void initializeParticleRandom();
};