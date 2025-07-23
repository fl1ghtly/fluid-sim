#include "FluidParameters.h"

FluidParameters::FluidParameters(
		float gravity, 
		float damping, 
		float restDensity, 
		float stiffness,
		float viscosity, 
		float smoothingRadius
	)
	: gravity(gravity),
	damping(damping),
	restDensity(restDensity),
	stiffness(stiffness),
	viscosity(viscosity),
	smoothingRadius(smoothingRadius),
	sqSmoothingRadius(smoothingRadius * smoothingRadius) {}