#include "Boundary.h"

Boundary::Boundary(float smoothingRadius) : smoothingRadius(smoothingRadius) {
	boundaryMass = INFINITY;
	isStatic = true;
	boundaryVel = {0.f, 0.f};
	spatialGrid = createHashmap<GridCell, std::vector<std::size_t>>();
}

Boundary::Boundary(
	float smoothingRadius, 
	float mass, 
	Vector2f initialVel
) : smoothingRadius(smoothingRadius),
boundaryMass(mass),
boundaryVel(initialVel) {
	isStatic = false;
	spatialGrid = createHashmap<GridCell, std::vector<std::size_t>>();
}

void Boundary::applyForceAndTorque(Vector2f force, float dt) {
	if (isStatic) return;
	boundaryVel += dt * force / boundaryMass;

	synchronizeBoundaryParticles(dt);
}

std::size_t Boundary::getNumBoundaryParticles() const {
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

void Boundary::createPolygon(std::vector<Vector2f> vertices, float compression) {
	const float particleSpacing = 1.f / std::clamp(compression, 1E-6f, 1.f);
	centerOfMass = {0.f, 0.f};
	// Construct lines between each vertex in sequence
	// The last vertex connects to the first vertex in the sequence
	for (std::size_t i = 0; i < vertices.size(); i++) {
		const Vector2f start = vertices[i];
		const Vector2f end = vertices[(i + 1) % vertices.size()];
		const Vector2f line = end - start;
		
		const float lineParticleCount = std::floor(line.magnitude() / particleSpacing) + 1.f;
		
		const Vector2f norm = line.normalize();
		const float dx = norm.x;
		const float dy = norm.y;
		
		for (float j = 0.f; j < lineParticleCount; j++) {
			particlePos.push_back({
				start.x + particleSpacing * j * dx,
				start.y + particleSpacing * j * dy 
			});
		}

		centerOfMass += start;
	}
	centerOfMass /= static_cast<float>(vertices.size());
}

void Boundary::createBox(Vector2f corner1, Vector2f corner2, float compression) {
	const float particleSpacing = 1.f / compression;
	// Only care about magnitude of width/height
	const float width = std::fabs(corner2.x - corner1.x);
	const float height = std::fabs(corner2.y - corner1.y);
	
	// Get amount of particles in the x and y direction
	const float widthParticles = std::floor(width / particleSpacing) + 1.f;
	const float heightParticles = std::floor(height / particleSpacing) + 1.f;

	// Calculate direction of line depending the relative position of the corners
	const float widthDir = (corner1.x < corner2.x) ? 1.f : -1.f;
	const float heightDir = (corner1.y < corner2.y) ? 1.f : -1.f;

	// Construct top and bottom of box
	for (float i = 0.f; i < widthParticles; i++) {
		const Vector2f top = {
			corner1.x + particleSpacing * i  * widthDir, 
			corner1.y
		};
		const Vector2f bottom = {
			corner1.x + particleSpacing * i * widthDir, 
			corner2.y
		};
		particlePos.push_back(top);
		particlePos.push_back(bottom);
	}

	// Construct left and right of box
	for (float i = 0.f; i < heightParticles; i++) {
		const Vector2f left = {
			corner1.x, 
			corner1.y + particleSpacing * i * heightDir
		};
		const Vector2f right = {
			corner2.x, 
			corner1.y + particleSpacing * i * heightDir
		};
		particlePos.push_back(left);
		particlePos.push_back(right);
	}
	centerOfMass = {corner1.x + width / 2.f, corner1.y + height / 2.f};
}

void Boundary::createCircle(Vector2f origin, float radius, float compression) {
	const float particleSpacing = 1.f / compression;
	constexpr float PI = std::numbers::pi_v<float>;
	const float circleParticles = std::floor(2.f * PI * radius / particleSpacing);

	for (float i = 0; i < circleParticles; i++) {
		const float deg = 360.f * i / circleParticles;
		const float radian = deg * PI / 180.f;
		particlePos.push_back({
			origin.x + radius * cosf(radian),
			origin.y + radius * sinf(radian)
		});
	}
	centerOfMass = origin;
}

void Boundary::activateBoundary() {
	calculateParticleVolume();
}

void Boundary::calculateParticleVolume() {
	findNeighbors();
	particleVolume.resize(particlePos.size());
	for (std::size_t i = 0; i < particlePos.size(); i++) {
		float sum = 0.f;
		for (std::size_t j : neighbors[i]) {
			const float dist = (particlePos[i] - particlePos[j]).magnitude();
			sum += poly6Kernel(dist);
		}
		particleVolume[i] = 1 / sum;
	}
}

void Boundary::synchronizeBoundaryParticles(float dt) {
	for (std::size_t i = 0; i < particlePos.size(); i++) {
		particlePos[i] += dt * boundaryVel;
	}
}

void Boundary::buildSpatialGrid() {
	// Cell Size = Kernel Support Length is optimal
	const float cellSize = smoothingRadius;

	spatialGrid->clear();
	#pragma omp parallel for
	for (std::size_t i = 0; i < particlePos.size(); i++) {
		GridCell c = {particlePos[i], cellSize};
		if (spatialGrid->find(c)) {
			#pragma omp critical
			spatialGrid->get(c).push_back(i);
		} else {
			spatialGrid->insert(c, {i});
		}
	}
}

void Boundary::findNeighbors() {
	const float cellSize = smoothingRadius;
	buildSpatialGrid();
	neighbors.resize(particlePos.size());

	#pragma omp parallel for
	for (std::size_t i = 0; i < particlePos.size(); i++) {
		std::vector<std::size_t> neighborIndices;
		const GridCell center = {particlePos[i], cellSize};
		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				const GridCell cell = {
					static_cast<int>(center.x) + dx, 
					static_cast<int>(center.y) + dy, 
					cellSize
				};
				const std::vector<std::size_t> cellIndices = spatialGrid->get(cell);

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
	if (dist < 0.f || dist > smoothingRadius) return 0.f;
	constexpr float PI = std::numbers::pi_v<float>;
	const float factor = smoothingRadius * smoothingRadius - dist * dist;
	const float poly6C = 4.f / (PI * pow(smoothingRadius, 8.f));
	return poly6C * factor * factor * factor;
}