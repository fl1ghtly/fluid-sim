#pragma once
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <unordered_map>
#include <random>
#include "Vector2f.h"
#include "GridCell.h"

class Fluid {
	private:
		int width;
		int height;
		int numParticles;
		float radius;
		float fluidDensity;
		std::vector<Vector2f> position;
		std::vector<Vector2f> velocity;
		std::vector<float> mass;
		std::vector<float> density;
		std::vector<float> pressure;
		float smoothingLen;
		std::unordered_map<GridCell, std::vector<int>> spatialGrid;
		std::vector<std::vector<int>> neighbors;

		void initializeParticleValues();
		void calculateDensity(float dt);
		void calculatePressure();
		std::vector<float> calculateAlphaFactors();
		void applyNonPressureForce(float dt);
		void applyPressureForce(float dt);
		void applyBoundaryCondition();
		void correctDensityError(std::vector<float> alpha, float dt);
		void correctDivergenceError(std::vector<float> alpha, float dt);
		float calculateTimeStep();
		void buildSpatialGrid();
		void findNeighbors();
		float smoothingKernel(float dist, float smoothingLength);
		Vector2f smoothingGradient(Vector2f r, float smoothingLength);
		template <typename T>
		Vector2f gradient(int particleIndex, std::vector<T> field);
		float divergence(int particleIndex, std::vector<Vector2f> field);
		Vector2f laplacian(int particleIndex, std::vector<Vector2f> field);
	public:
		Fluid(int width, int heigh, int numParticles, float particleRadius, float density);
		void update();
		std::vector<Vector2f> getPosition();
		std::vector<Vector2f> getVelocity();
		void initializeParticleGrid(int x, int y, int width, int height);
		void initializeParticleRandom();
};