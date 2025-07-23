#include "FluidParameters.h"

FluidParameters::FluidParameters(
		Vector2f downDirection, 
		float gravity, 
		float damping, 
		float restDensity, 
		float stiffness,
		float viscosity, 
		float smoothingRadius
	)
	: downDirection(downDirection),
	gravity(gravity),
	damping(damping),
	restDensity(restDensity),
	stiffness(stiffness),
	viscosity(viscosity),
	smoothingRadius(smoothingRadius),
	sqSmoothingRadius(smoothingRadius * smoothingRadius) {}

Vector2f FluidParameters::getGravityVector() const {
	return downDirection * gravity;
}