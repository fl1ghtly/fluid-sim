#pragma once
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <random>
#include <ranges>
#include "Vector2f.h"
#include "GridCell.h"
#include "FluidParameters.h"
#include "Boundary.h"

#ifdef TRACY_ENABLED
#include <tracy/Tracy.hpp>
#define SimZoneScoped ZoneScoped
#else
#define SimZoneScoped
#endif

#ifdef MULTIPROCESSING_ENABLED
#include "Hashmap/HashmapMP.h"
#else
#include "Hashmap/HashmapST.h"
#endif

class Simulation {
	public:
		Simulation(
			float width, 
			float heigh, 
			FluidParameters& params, 
			Vector2f downDir={0.f, 1.f}, 
			float fixedTimestep=1E-2f
		);
		void update();
		std::vector<Vector2f> getPosition();
		std::vector<Vector2f> getVelocity();
		std::vector<float> getPressure();
		std::vector<float> getDensity();
		std::size_t getNumParticles();
		void setDownDirection(Vector2f d);
		float getPressureAtPoint(const Vector2f pos);
		void initializeParticleGrid(Vector2f center, std::size_t gridWidth, std::size_t amount);
		void initializeParticleRandom(Vector2f min, Vector2f max, std::size_t amount);
		void initializeParticleCircle(Vector2f origin, float radius);
		void addBoundary(Boundary b);
		void addBoundary(std::vector<Boundary> b);

	private:
		void initializeParticleValues(std::size_t amount);
		void calculateDensity(float dt);
		void calculatePressure();
		void applyNonPressureForce(float dt);
		void applyPressureForce(float dt);
		void applyForceToBoundary(float dt);
		void applyWorldBoundary();
		float calculateTimeStep();
		void buildSpatialGrid();
		void sortZIndex();
		void findNeighbors();
		float poly6Kernel(float dist, float smoothingLength);
		Vector2f poly6Gradient(Vector2f r, float smoothingLength);
		Vector2f spikyGradient(Vector2f r, float smoothingLength);
		float viscosityLaplacian(float dist, float smoothingLength);

		float width;
		float height;
		int currentStep;
		std::size_t numParticles;
		FluidParameters& params;
		Vector2f downDir;
		float fixedTimestep;
		std::vector<Vector2f> position;
		std::vector<Vector2f> velocity;
		std::vector<float> mass;
		std::vector<float> density;
		std::vector<float> pressure;
		Hashmap<GridCell, std::vector<std::size_t>>* spatialGrid;
		std::vector<std::vector<std::size_t>> neighbors;
		std::vector<Boundary> boundaries;
		std::vector<Vector2f> boundaryPos;
		std::vector<float> boundaryVolume;
		std::vector<std::size_t> boundaryStartIndex;
		std::vector<Vector2f> boundaryForce;
		Hashmap<GridCell, std::vector<std::size_t>>* boundarySpatialGrid;
		std::vector<std::vector<std::size_t>> boundaryNeighbors;
		float poly6C, spikyGC, viscosityLC;
};