#include "fluid.h"

#define GRAVITY 9.8f
#define DAMPING 0.95f

constexpr Vector2f VecDown = Vector2f(0.f, 1.f);

Fluid::Fluid(int width, int height, int numParticles, float smoothingRadius, float density) 
	: width(width), 
	  height(height), 
	  numParticles(numParticles),
	  smoothingLen(smoothingRadius), 
	  fluidDensity(density), 
	  position(numParticles), 
	  velocity(numParticles), 
	  mass(numParticles), 
	  density(numParticles), 
	  pressure(numParticles),
	  neighbors(numParticles)
{
}

void Fluid::initializeParticleValues() {
	const float particleArea = (4.f / 9.f) * smoothingLen * smoothingLen;
	const float particleMass = fluidDensity * particleArea;
	std::fill(velocity.begin(), velocity.end(), Vector2f(0.f, 0.f));
	std::fill(mass.begin(), mass.end(), particleMass);
	std::fill(density.begin(), density.end(), fluidDensity);
	std::fill(pressure.begin(), pressure.end(), 0.f);
}

void Fluid::initializeParticleGrid(int gridWidth) {
	constexpr float spacingFactor = 2.f;
	const float spacing = smoothingLen / spacingFactor;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> jitter(0, spacing);
	
	Vector2f center(width / 2, height / 2);
	center.x -= 2 * gridWidth;
	center.y -= 2 * numParticles / gridWidth;
	for (int i = 0; i < numParticles; i++) {
		position[i] = Vector2f(center.x + jitter(gen) + spacing * (i % gridWidth), 
							   center.y + jitter(gen) + spacing * (i / gridWidth));
	}
	initializeParticleValues();
}

void Fluid::initializeParticleRandom() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> xDistr(0, width);
	std::uniform_real_distribution<float> yDistr(0, height);
	for (int i = 0; i < numParticles; i++) {
		position[i] = Vector2f(xDistr(gen), yDistr(gen));
	}
	initializeParticleValues();
}

void Fluid::update() {
	const float deltaTime = 0.01f;
	// const float deltaTime = calculateTimeStep();
	findNeighbors();
	applyNonPressureForce(deltaTime);
	calculateDensity(deltaTime);
	calculatePressure();
	applyPressureForce(deltaTime);       
	applyBoundaryCondition();
}

void Fluid::calculateDensity(float dt) {
	for (int i = 0; i < numParticles; i++) {
		density[i] = 0.f;
		for (auto j : neighbors[i]) {
			const Vector2f dist = position[i] - position[j];
			const float influence = poly6Kernel(dist.magnitude(), smoothingLen);
			const Vector2f velDiff = velocity[i] - velocity[j];
			const Vector2f smoothGrad = poly6Gradient(dist, smoothingLen);
			density[i] += mass[j] * influence + dt * velDiff.dot(smoothGrad);
		}
	}
}

void Fluid::calculatePressure() {
	constexpr float stiffness = 10000.f;
	for (int i = 0; i < numParticles; i++) {
		pressure[i] = stiffness * (pow(density[i] / fluidDensity, 7.f) - 1);
	}
}

void Fluid::applyNonPressureForce(float dt) {
	const float viscosity = 1E+6f;
	for (int i = 0; i < numParticles; i++) {
		// Force due to gravity
		Vector2f gForce(0.f, 0.f); 
		gForce = VecDown * GRAVITY * mass[i];
		
		// Force due to viscosity
		Vector2f vForce(0.f, 0.f);
		for (int j : neighbors[i]) {
			const float volume = mass[j] / density[j];
			const Vector2f velDiff = velocity[j] - velocity[i];
			const float vLaplacian = viscosityLaplacian((position[i] - position[j]).magnitude(), smoothingLen);
			vForce += volume * velDiff * vLaplacian;
		}
		vForce *= viscosity;
		
		const Vector2f Force = gForce + vForce;
		velocity[i] += dt * Force / mass[i];
		position[i] += dt * velocity[i];
	}
}

void Fluid::applyPressureForce(float dt) {
	std::vector<Vector2f> forces(numParticles);
	for (int i = 0; i < numParticles; i++) {
		Vector2f force(0.f, 0.f);
		for (int j : neighbors[i]) {
			const float volume = mass[j] / density[j];
			const float avgPressure = (pressure[i] + pressure[j]) / 2.f;
			const Vector2f grad = spikyGradient(position[i] - position[j], smoothingLen);
			force += -volume * avgPressure * grad;
		}
		forces[i] =  force;
	}
	
	for (int i = 0; i < numParticles; i++) {
		velocity[i] += dt * forces[i] / mass[i];
		position[i] += dt * velocity[i];
	}
}

void Fluid::applyBoundaryCondition() {
	for (int i = 0; i < numParticles; i++) {
		auto& p = position[i];
		auto& v = velocity[i];

		if (p.x <= 0.f) {
			p.x = 0.f;
			v.x *= -DAMPING;
		} else if (p.x >= width) {
			p.x = width - 1;
			v.x *= -DAMPING;
		}
		
		if (p.y <= 0.f) {
			p.y = 0.f;
			v.y *= -DAMPING;
		} else if (p.y >= height) {
			p.y = height - 1;
			v.y *= -DAMPING;
		}
	}
}

float Fluid::calculateTimeStep() {
	const Vector2f v = *std::max_element(velocity.begin(), velocity.end());
	const float vMax = v.magnitude();
	const float timestep = 0.4f * 2.f * smoothingLen / (vMax + 1E-6f);

	return std::clamp(timestep, 1E-6f, 1E-2f);
}

void Fluid::buildSpatialGrid() {
	// Cell Size = Kernel Support Length is optimal
	const float cellSize = smoothingLen;

	spatialGrid.clear();
	#pragma omp parallel for
	for (int i = 0; i < numParticles; i++) {
		GridCell c = {position[i], cellSize};
		spatialGrid[c].push_back(i);
	}
}

void Fluid::findNeighbors() {
	const float cellSize = smoothingLen;
	buildSpatialGrid();
	
	for (int i = 0; i < numParticles; i++) {
		std::vector<int> neighborIndices;
		const GridCell center = {position[i], cellSize};
		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				GridCell cell = {center.x + dx, center.y + dy, cellSize};
				const std::vector<int> cellIndices = spatialGrid[cell];
				for (auto neighbor : cellIndices) {
					// if (neighbor == i) continue;
					const float dist = (position[i] - position[neighbor]).magnitude();
					if (dist > smoothingLen) continue;
					neighborIndices.push_back(neighbor);
				}
			}
		}
		neighbors[i] = neighborIndices;
	}
}

template <typename T>
Vector2f Fluid::gradient(int particleIndex, std::vector<T> field) {
	Vector2f sum = Vector2f(0.f, 0.f);
	
	for (auto j : neighbors[particleIndex]) {
		if (particleIndex == j) continue;
		const float Ai = field[particleIndex] / (density[particleIndex] * density[particleIndex]);
		const float Aj = field[j] / (density[j] * density[j]);
		Vector2f grad = poly6Gradient(position[particleIndex] - position[j], smoothingLen);
		sum += mass[j] * (Ai + Aj) * grad;
		if (std::isnan(sum.x) || std::isnan(sum.y)) {
			printf("Ai %f, Aj %f, grad (%f, %f)\n", Ai, Aj, grad.x, grad.y);
		}
	}
	return density[particleIndex] * sum;
}

float Fluid::divergence(int particleIndex, std::vector<Vector2f> field) {
	float sum = 0.f;

	for (auto j : neighbors[particleIndex]) {
		if (particleIndex == j) continue;
		const Vector2f Aij = field[particleIndex] - field[j];
		const Vector2f grad = poly6Gradient(position[particleIndex] - position[j], smoothingLen);
		sum += mass[j] * Aij.dot(grad);
	}
	return -1.f / density[particleIndex] * sum;
}

Vector2f Fluid::laplacian(int particleIndex, std::vector<Vector2f> field) {
	Vector2f sum = Vector2f(0.f, 0.f);

	for (auto j : neighbors[particleIndex]) {
		if (particleIndex == j) continue;
		const float volume = mass[j] / density[j];
		const Vector2f Aij = field[particleIndex] - field[j];
		const Vector2f Xij = position[particleIndex] - position[j];
		const Vector2f smoothGrad = poly6Gradient(Xij, smoothingLen);
		const float numerator = Xij.dot(smoothGrad);
		const float denominator = Xij.dot(Xij) + 0.01 * smoothingLen * smoothingLen;
		sum += volume * Aij * numerator / denominator;
	}
	return 2.f * sum;
}

std::vector<Vector2f> Fluid::getPosition() {
	return position;
}

std::vector<Vector2f> Fluid::getVelocity() {
	return velocity;
}

float Fluid::getPressureAtPoint(const Vector2f point) {
	const float cellSize = smoothingLen;
	const GridCell center = {point, cellSize};
	float p = 0.f;

	for (int dx = -1; dx <= 1; dx++) {
		for (int dy = -1; dy <= 1; dy++) {
			GridCell cell = {center.x + dx, center.y + dy, cellSize};
			const std::vector<int> cellIndices = spatialGrid[cell];
			for (auto neighbor : cellIndices) {
				p += pressure[neighbor];			
			}
		}
	}
	return p;
}

float Fluid::poly6Kernel(float dist, float smoothingLength) {
	if (dist < 0 || dist > smoothingLength) return 0.f;
	const float coeff = 4 / (M_PI * pow(smoothingLength, 8.f));
	const float factor = smoothingLength * smoothingLength - dist * dist;
	return coeff * factor * factor * factor;
}

Vector2f Fluid::poly6Gradient(Vector2f r, float smoothingLength) {
	const float dist = r.magnitude();
	if (dist <= 0 || dist > smoothingLength) return Vector2f(0.f, 0.f);
	const float coeff = 4 / (M_PI * pow(smoothingLength, 8.f));
	const float factor = -6 * (smoothingLength * smoothingLength - dist * dist);
	return coeff * factor * factor * r;
}

Vector2f Fluid::spikyGradient(Vector2f r, float smoothingLength) {
	const float dist = r.magnitude();
	if (dist <= 0 || dist > smoothingLength) return Vector2f(0.f, 0.f);
	const float coeff = -30.f / (M_PI * pow(smoothingLength, 5.f));
	const float factor = smoothingLength - dist;
	return coeff * factor * factor * r / dist;
}

float Fluid::viscosityLaplacian(float dist, float smoothingLength) {
	if (dist <= 0 || dist > smoothingLength) return 0.f;
	const float coeff = 40.f / (M_PI * pow(smoothingLen, 5.f));
	const float factor = smoothingLength - dist;
	return coeff * factor;
}