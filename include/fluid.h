#pragma once
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <unordered_map>
#include <random>
#include "Vector2f.h"
#include "GridCell.h"
#include "FluidParameters.h"

class Fluid {
	private:
		int width;
		int height;
		FluidParameters& params;
		float fixedTimestep;
		std::vector<Vector2f> position;
		std::vector<Vector2f> velocity;
		std::vector<float> mass;
		std::vector<float> density;
		std::vector<float> pressure;
		std::unordered_map<GridCell, std::vector<int>> spatialGrid;
		std::vector<std::vector<int>> neighbors;

		void initializeParticleValues();
		void calculateDensity(float dt);
		void calculatePressure();
		void applyNonPressureForce(float dt);
		void applyPressureForce(float dt);
		void applyBoundaryCondition();
		float calculateTimeStep();
		void buildSpatialGrid();
		void findNeighbors();
		float poly6Kernel(float dist, float smoothingLength);
		Vector2f poly6Gradient(Vector2f r, float smoothingLength);
		Vector2f spikyGradient(Vector2f r, float smoothingLength);
		float viscosityLaplacian(float dist, float smoothingLength);
		template <typename T>
		Vector2f gradient(int particleIndex, std::vector<T> field);
		float divergence(int particleIndex, std::vector<Vector2f> field);
		Vector2f laplacian(int particleIndex, std::vector<Vector2f> field);
	public:
		Fluid(int width, int heigh, FluidParameters& params, float fixedTimestep=1E-2f);
		void update();
		std::vector<Vector2f> getPosition();
		std::vector<Vector2f> getVelocity();
		float getPressureAtPoint(const Vector2f pos);
		void initializeParticleGrid(int gridWidth);
		void initializeParticleRandom();
};