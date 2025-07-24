#pragma once
#include <vector>
#include <oneapi/tbb/concurrent_hash_map.h>
#include "Vector2f.h"
#include "GridCell.h"

class Boundary {
	public:
		Boundary(float smoothingRadius);
		Boundary(
			float smoothingRadius,
			float mass,
			Vector2f initialVel={0.f, 0.f}
		);
		void createPolygon(std::vector<Vector2f> vertices);
		void createBox(Vector2f topLeft, Vector2f bottomRight);
		void createCircle(Vector2f origin, float radius);
		void applyForce(Vector2f force, float dt);
		int getNumBoundaryParticles() const;
		std::vector<Vector2f> getBoundaryParticlePositions() const;
		std::vector<float> getBoundaryParticleVolume() const;
		Vector2f getCenterPosition() const;
		
	private:
		float boundaryMass;
		float smoothingRadius;
		bool isStatic;
		Vector2f boundaryVel;
		float particleSpacing;
		Vector2f centerOfMass;
		std::vector<Vector2f> particlePos;
		std::vector<float> particleVolume;
		tbb::concurrent_hash_map<GridCell, std::vector<int>> spatialGrid;
		std::vector<std::vector<int>> neighbors;

		void calculateParticleVolume();
		void synchronizeBoundaryParticles(float dt);
		void buildSpatialGrid();
		void findNeighbors();
		float poly6Kernel(float dist);
};