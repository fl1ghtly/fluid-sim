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
		float sqSmoothingRadius;
		
		FluidParameters(
			Vector2f downDirection, 
			float gravity, 
			float damping, 
			float restDensity, 
			float stiffness, 
			float viscosity, 
			float smoothingRadius
		);

		Vector2f getGravityVector() const;
};