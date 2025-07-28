#pragma once
#include <vector>
#include "Vector2f.h"
#include "GridCell.h"

#ifdef MULTIPROCESSING_ENABLED
#include "Hashmap/HashmapMP.h"
#else
#include "Hashmap/HashmapST.h"
#endif

class Boundary {
	public:
		Boundary(float smoothingRadius);
		Boundary(
			float smoothingRadius,
			float mass,
			Vector2f initialVel={0.f, 0.f}
		);
		void createPolygon(std::vector<Vector2f> vertice, float compression=1.f);
		void createBox(Vector2f topLeft, Vector2f bottomRight, float compression=1.f);
		void createCircle(Vector2f origin, float radius, float compression=1.f);
		void activateBoundary();
		void applyForceAndTorque(Vector2f force, float dt);
		int getNumBoundaryParticles() const;
		std::vector<Vector2f> getBoundaryParticlePositions() const;
		std::vector<float> getBoundaryParticleVolume() const;
		Vector2f getCenterPosition() const;
		
	private:
		float boundaryMass;
		float smoothingRadius;
		bool isStatic;
		Vector2f boundaryVel;
		Vector2f centerOfMass;
		std::vector<Vector2f> particlePos;
		std::vector<float> particleVolume;
		Hashmap<GridCell, std::vector<int>>* spatialGrid;
		std::vector<std::vector<int>> neighbors;

		void calculateParticleVolume();
		void synchronizeBoundaryParticles(float dt);
		void buildSpatialGrid();
		void findNeighbors();
		float poly6Kernel(float dist);
};