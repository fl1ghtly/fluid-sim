#pragma once
#include "Vector2f.h"

class FluidParameters {
	public:
		float gravity; 
		float damping;
		float restDensity;
		float stiffness;
		float viscosity;
		float smoothingRadius;
		float sqSmoothingRadius;
		
		FluidParameters(
			float gravity, 
			float damping, 
			float restDensity, 
			float stiffness, 
			float viscosity, 
			float smoothingRadius
		);
};