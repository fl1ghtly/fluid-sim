#include "Boundary.h"

constexpr float spacingFactor = 2.5f;

Boundary::Boundary(
	std::vector<Vector2f> vertices, 
	float mass, 
	float smoothingRadius,
	Vector2f initialVel
	) : boundaryMass(mass),
	smoothingRadius(smoothingRadius),
	boundaryVel(initialVel)
	{
		particleSpacing = smoothingRadius / spacingFactor;
		initializeBoundaryParticles(vertices);
	}

Boundary::Boundary(
	Vector2f topLeft, 
	Vector2f bottomRight, 
	float mass, 
	float smoothingRadius,
	Vector2f initialVel
	) : boundaryMass(mass), 
	smoothingRadius(smoothingRadius),
	boundaryVel(initialVel)
	{
		particleSpacing = smoothingRadius / spacingFactor;
		initializeBoundaryParticles(topLeft, bottomRight);
	}

Boundary::Boundary(
	Vector2f origin, 
	float radius, 
	float mass, 
	float smoothingRadius,
	Vector2f initialVel
	) : boundaryMass(mass),
	smoothingRadius(smoothingRadius),
	boundaryVel(initialVel)
	{
		particleSpacing = smoothingRadius / spacingFactor;
		initializeBoundaryParticles(origin, radius);
	}

void Boundary::applyForce(Vector2f force, float dt) {
	boundaryVel += dt * force / boundaryMass;

	synchronizeBoundaryParticles(dt);
}

int Boundary::getNumBoundaryParticles() const {
	return particlePos.size();
}

std::vector<Vector2f> Boundary::getBoundaryParticlePositions() const {
	return particlePos;
}

Vector2f Boundary::getCenterPosition() const {
	return centerOfMass;
}

void Boundary::initializeBoundaryParticles(std::vector<Vector2f> vertices) {
	// Polygon

	centerOfMass = {0.f, 0.f};
	// Construct lines between each vertex in sequence
	// The last vertex connects to the first vertex in the sequence
	for (int i = 0; i < vertices.size(); i++) {
		const Vector2f start = vertices[i];
		const Vector2f end = vertices[(i + 1) % vertices.size()];
		const Vector2f line = end - start;
		
		const int lineParticleCount = (int) line.magnitude() / particleSpacing + 1;
		
		const Vector2f norm = line.normalize();
		const float dx = norm.x;
		const float dy = norm.y;
		
		for (int i = 0; i < lineParticleCount; i++) {
			particlePos.push_back({
				start.x + particleSpacing * i * dx,
				start.y + particleSpacing * i * dy 
			});
		}

		centerOfMass += start;
	}

	centerOfMass /= vertices.size();
}

void Boundary::initializeBoundaryParticles(Vector2f topLeft, Vector2f bottomRight) {
	// Box 
	
	// Only care about magnitude of width/height
	const float width = bottomRight.x - topLeft.x;
	const float height = bottomRight.y - topLeft.y;
	
	// Get amount of particles in the x and y direction
	const int widthParticles = (int) width / particleSpacing + 1;
	const int heightParticles = (int) height / particleSpacing + 1;
	
	// Construct top and bottom of box
	for (int i = 0; i < widthParticles; i++) {
		const Vector2f top = {topLeft.x + particleSpacing * i, topLeft.y};
		const Vector2f bottom = {topLeft.x + particleSpacing * i, bottomRight.y};
		particlePos.push_back(top);
		particlePos.push_back(bottom);
	}

	// Construct left and right of box
	for (int i = 0; i < heightParticles; i++) {
		const Vector2f left = {topLeft.x, topLeft.y + particleSpacing * i};
		const Vector2f right = {bottomRight.x, topLeft.y + particleSpacing * i};
		particlePos.push_back(left);
		particlePos.push_back(right);
	}

	centerOfMass = {topLeft.x + width / 2.f, topLeft.y + height / 2.f};
}

void Boundary::initializeBoundaryParticles(Vector2f origin, float radius) {
	// Circle
	
	const int circleParticles = (int) 2.f * M_PI * radius / particleSpacing;

	for (int i = 0; i < circleParticles; i++) {
		const float deg = 360.f * i / circleParticles;
		const float radian = deg * M_PI / 180.f;
		particlePos.push_back({
			origin.x + radius * cosf(radian),
			origin.y + radius * sinf(radian)
		});
	}

	centerOfMass = origin;
}

void Boundary::synchronizeBoundaryParticles(float dt) {
	for (int i = 0; i < particlePos.size(); i++) {
		particlePos[i] += dt * boundaryVel;
	}
}