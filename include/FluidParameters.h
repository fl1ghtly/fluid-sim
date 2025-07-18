#pragma once
#include "Vector2f.h"

class FluidParameters {
	public:
		Vector2f downDirection;
		float gravity; 
		float damping;
		float restDensity;
		float stiffness;
		float viscosity;
		float smoothingRadius;
		int numParticles;
		
		FluidParameters(Vector2f downDirection, float gravity, float damping, 
						float restDensity, float stiffness, float viscosity, 
						float smoothingRadius, int numParticles);

		Vector2f getGravityVector() const;
};