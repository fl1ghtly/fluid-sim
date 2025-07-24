#include "Boundary.h"

constexpr float spacingFactor = 2.5f;

Boundary::Boundary(float smoothingRadius) : smoothingRadius(smoothingRadius) {
	boundaryMass = INFINITY;
	isStatic = true;
	boundaryVel = {0.f, 0.f};
	particleSpacing = smoothingRadius / spacingFactor;
}

Boundary::Boundary(
	float smoothingRadius, 
	float mass, 
	Vector2f initialVel
) : smoothingRadius(smoothingRadius),
boundaryMass(mass),
boundaryVel(initialVel) {
	isStatic = false;
	particleSpacing = smoothingRadius / spacingFactor;
}

void Boundary::applyForce(Vector2f force, float dt) {
	if (isStatic) return;
	boundaryVel += dt * force / boundaryMass;

	synchronizeBoundaryParticles(dt);
}

int Boundary::getNumBoundaryParticles() const {
	return particlePos.size();
}

std::vector<Vector2f> Boundary::getBoundaryParticlePositions() const {
	return particlePos;
}

std::vector<float> Boundary::getBoundaryParticleVolume() const {
	return particleVolume;
}

Vector2f Boundary::getCenterPosition() const {
	return centerOfMass;
}

void Boundary::createPolygon(std::vector<Vector2f> vertices) {
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

	calculateParticleVolume();
}

void Boundary::createBox(Vector2f topLeft, Vector2f bottomRight) {
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

	calculateParticleVolume();
}

void Boundary::createCircle(Vector2f origin, float radius) {
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

	calculateParticleVolume();
}

void Boundary::calculateParticleVolume() {
	findNeighbors();
	particleVolume.resize(particlePos.size());
	for (int i = 0; i < particlePos.size(); i++) {
		float sum = 0.f;
		for (int j : neighbors[i]) {
			const float dist = (particlePos[i] - particlePos[j]).magnitude();
			sum += poly6Kernel(dist);
		}
		particleVolume[i] = 1 / sum;
	}
}

void Boundary::synchronizeBoundaryParticles(float dt) {
	for (int i = 0; i < particlePos.size(); i++) {
		particlePos[i] += dt * boundaryVel;
	}
}

void Boundary::buildSpatialGrid() {
	// Cell Size = Kernel Support Length is optimal
	const float cellSize = smoothingRadius;

	spatialGrid.clear();
	#pragma omp parallel for
	for (int i = 0; i < particlePos.size(); i++) {
		GridCell c = {particlePos[i], cellSize};
		tbb::concurrent_hash_map<GridCell, std::vector<int>>::accessor a;
		if (spatialGrid.insert(a, c)) {
			a->second = {i};
		} else {
			a->second.push_back(i);
		}
		a.release();
	}
}

void Boundary::findNeighbors() {
	const float cellSize = smoothingRadius;
	buildSpatialGrid();
	neighbors.resize(particlePos.size());

	#pragma omp parallel for
	for (int i = 0; i < particlePos.size(); i++) {
		std::vector<int> neighborIndices;
		const GridCell center = {particlePos[i], cellSize};
		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				const GridCell cell = {center.x + dx, center.y + dy, cellSize};
				std::vector<int> cellIndices;
				{
					tbb::concurrent_hash_map<GridCell, std::vector<int>>::const_accessor ca;
					if (spatialGrid.find(ca, cell)) cellIndices = ca->second;
				}

				for (auto neighbor : cellIndices) {
					const Vector2f rij = particlePos[i] - particlePos[neighbor];
					const float dist = rij.magnitude();
					if (dist > smoothingRadius) continue;
					neighborIndices.push_back(neighbor);
				}
			}
		}
		neighbors[i] = std::move(neighborIndices);
	}
}

float Boundary::poly6Kernel(float dist) {
	if (dist < 0 || dist > smoothingRadius) return 0.f;
	const float factor = smoothingRadius * smoothingRadius - dist * dist;
	const float poly6C = 4 / (M_PI * pow(smoothingRadius, 8.f));
	return poly6C * factor * factor * factor;
}