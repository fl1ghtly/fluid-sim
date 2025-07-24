#pragma once
#include <vector>
#include "Vector2f.h"

class Boundary {
	public:
		Boundary(
			std::vector<Vector2f> vertices, 
			float mass, 
			float smoothingRadius,
			Vector2f initialVel={0.f, 0.f} 
		);
		Boundary(
			Vector2f topLeft, 
			Vector2f bottomRight, 
			float mass, 
			float smoothingRadius,
			Vector2f initialVel={0.f, 0.f} 
		);
		Boundary(
			Vector2f origin, 
			float radius, 
			float mass, 
			float smoothingRadius,
			Vector2f initialVel={0.f, 0.f} 
		);
		void applyForce(Vector2f force, float dt);
		int getNumBoundaryParticles() const;
		std::vector<Vector2f> getBoundaryParticlePositions() const;
		Vector2f getCenterPosition() const;
	private:
		float boundaryMass;
		float smoothingRadius;
		Vector2f boundaryVel;
		float particleSpacing;
		Vector2f centerOfMass;
		std::vector<Vector2f> particlePos;

		void initializeBoundaryParticles(std::vector<Vector2f> vertices);
		void initializeBoundaryParticles(Vector2f topLeft, Vector2f bottomRight);
		void initializeBoundaryParticles(Vector2f origin, float radius);
		void synchronizeBoundaryParticles(float dt);
};