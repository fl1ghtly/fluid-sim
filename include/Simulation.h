#pragma once
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <random>
#include <ranges>
#include <oneapi/tbb/concurrent_hash_map.h>
#include "Vector2f.h"
#include "GridCell.h"
#include "FluidParameters.h"

#include <tracy/Tracy.hpp>

class Simulation {
	private:
		int width;
		int height;
		int currentStep;
		FluidParameters& params;
		float fixedTimestep;
		std::vector<Vector2f> position;
		std::vector<Vector2f> velocity;
		std::vector<float> mass;
		std::vector<float> density;
		std::vector<float> pressure;
		tbb::concurrent_hash_map<GridCell, std::vector<int>> spatialGrid;
		std::vector<std::vector<int>> neighbors;
		float poly6C, spikyGC, viscosityLC;

		void initializeParticleValues();
		void calculateDensity(float dt);
		void calculatePressure();
		void applyNonPressureForce(float dt);
		void applyPressureForce(float dt);
		void applyBoundaryCondition();
		float calculateTimeStep();
		void buildSpatialGrid();
		void sortZIndex();
		void findNeighbors();
		float poly6Kernel(float dist, float smoothingLength);
		Vector2f poly6Gradient(Vector2f r, float smoothingLength);
		Vector2f spikyGradient(Vector2f r, float smoothingLength);
		float viscosityLaplacian(float dist, float smoothingLength);
	public:
		Simulation(int width, int heigh, FluidParameters& params, float fixedTimestep=1E-2f);
		void update();
		std::vector<Vector2f> getPosition();
		std::vector<Vector2f> getVelocity();
		std::vector<float> getPressure();
		std::vector<float> getDensity();
		float getPressureAtPoint(const Vector2f pos);
		void initializeParticleGrid(float x, float y, int gridWidth);
		void initializeParticleRandom();
};