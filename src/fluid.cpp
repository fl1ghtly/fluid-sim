#include "fluid.h"

#define GRAVITY 9.8
#define DAMPING 0.5f

constexpr Vector2f VecDown = Vector2f(0.f, 1.f);

Fluid::Fluid(int width, int height, int numParticles, float particleRadius, float density) 
	: width(width), 
	  height(height), 
	  numParticles(numParticles),
	  radius(particleRadius), 
	  fluidDensity(density), 
	  position(numParticles), 
	  velocity(numParticles), 
	  mass(numParticles), 
	  density(numParticles), 
	  pressure(numParticles),
	  neighbors(numParticles)
{
	smoothingLen = 2.f * radius;
}

void Fluid::initializeParticleValues() {
	const float particleArea = radius * radius;
	const float particleMass = fluidDensity * particleArea;
	std::fill(velocity.begin(), velocity.end(), Vector2f(0.f, 0.f));
	std::fill(mass.begin(), mass.end(), particleMass);
	std::fill(density.begin(), density.end(), fluidDensity);
	std::fill(pressure.begin(), pressure.end(), 0.f);
	// calculateDensity();
	// calculatePressure();
}

void Fluid::initializeParticleGrid(int x, int y, int w, int h) {
	const float spacingFactor = 1.5f;
	const float spacing = smoothingLen / spacingFactor;
	for (int i = 0; i < numParticles; i++) {
		position[i] = Vector2f(x + spacing * (i % w), y + spacing * (i / w));
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
	// float deltaTime = 0.4 * smoothingLen / (calculateVmax() + 1E-6);
	// float deltaTime = calculateTimeStep();
	float deltaTime = 0.1f;
	buildSpatialGrid();
	findNeighbors();
	applyNonPressureForce(deltaTime);
	calculateDensity(deltaTime);
	calculatePressure();
	applyPressureForce(deltaTime);       
	applyBoundaryCondition();
}

std::vector<float> Fluid::calculateAlphaFactors() {
	std::vector<float> alpha(numParticles);
	for (int i = 0; i < numParticles; i++) {
		Vector2f sumLeft = Vector2f(0.f, 0.f);
		float sumRight = 0.f;
		for (int j = 0; j < numParticles; j++) {
			if (i == j) continue;
			const Vector2f mWij = mass[j] * smoothingGradient(position[i] - position[j], smoothingLen);
			sumLeft += mWij;
			sumRight += mWij.magnitude() * mWij.magnitude();
		}
		const float lhs = sumLeft.magnitude() * sumLeft.magnitude();
		alpha[i] = std::max((float)10E-6, density[i] / (lhs + sumRight));
	}
	return alpha;
}

void Fluid::calculateDensity(float dt) {
	for (int i = 0; i < numParticles; i++) {
		density[i] = 0.f;
		// std::vector<int> neighbors = findNeighbors(i);
		for (auto j : neighbors[i]) {
			// if (i == j) continue;
			const Vector2f dist = position[i] - position[j];
			const float influence = smoothingKernel(dist.magnitude(), smoothingLen);
			const Vector2f velDiff = velocity[i] - velocity[j];
			const Vector2f smoothGrad = smoothingGradient(dist, smoothingLen);
			density[i] += mass[j] * influence + dt * velDiff.dot(smoothGrad);
		}
	}
}

void Fluid::calculatePressure() {
	constexpr float stiffness = 100.f;
	for (int i = 0; i < numParticles; i++) {
		// Rest density for water is 1000 kg/m^3
		pressure[i] = stiffness * (std::pow(density[i] / fluidDensity, 7.f) - 1);
	}
}

void Fluid::applyNonPressureForce(float dt) {
	const float viscosity = 10E-6;
	for (int i = 0; i < numParticles; i++) {
		Vector2f Force = {0.f, 0.f};
		// Force due to gravity
		// Force += VecDown * GRAVITY * mass[i];
		//Force due to viscosity
		const Vector2f lap = laplacian(i, velocity);
		Force += mass[i] * viscosity * laplacian(i, velocity);
		
		// if (std::isnan(Force.x) || std::isnan(Force.y)) {
		// 	printf("grad (%f, %f)\n", grad.x, grad.y);
		// }
		velocity[i] += dt * Force / mass[i];
		position[i] += dt * velocity[i];
	}
}

void Fluid::applyPressureForce(float dt) {
	for (int i = 0; i < numParticles; i++) {
		Vector2f Force = {0.f, 0.f};
		const Vector2f grad = gradient(i, pressure);
		Force += -mass[i] / density[i] * grad;
		
		velocity[i] += dt * Force / mass[i];
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

void Fluid::correctDivergenceError(std::vector<float> alpha, float dt) {
	std::vector<float> densityChange;
	// TODO loop until divergence error less than threshold
	for (int i = 0; i < numParticles; i++) {
		const float div = divergence(i, velocity);
		densityChange.push_back(-density[i] * div);
	}

	// TODO optimize
	for (int i = 0; i < numParticles; i++) {
		Vector2f sum = Vector2f(0.f, 0.f);
		const float stiffnessI = 1.f / dt * densityChange[i] * alpha[i];
		for (int j = 0; j < numParticles; j++) {
			if (i == j) continue;
			const float stiffnessJ = 1.f / dt * densityChange[j] * alpha[j];
			const float stiffIJ = stiffnessI / density[i] + stiffnessJ / density[j];
			// TODO double check
			sum += mass[j] * stiffIJ * smoothingGradient(position[i] - position[j], smoothingLen);
		}
		velocity[i] -= dt * sum;
	}
}

void Fluid::correctDensityError(std::vector<float> alpha, float dt) {
	// TODO loop until density error greater than some threshold
	// calculateDensity();
	for (int i = 0; i < numParticles; i++) {
		Vector2f sum = Vector2f(0.f, 0.f);
		const float stiffnessI = alpha[i] * (density[i] - fluidDensity) / (dt * dt);
		for (int j = 0; j < numParticles; j++) {
			if (i == j) continue;
			const float stiffnessJ = alpha[j] * (density[j] - fluidDensity) / (dt * dt);
			const float stiffIJ = stiffnessI / density[i] + stiffnessJ / density[j];
			sum += mass[j] * stiffIJ * smoothingGradient(position[i] - position[j], smoothingLen);
		}
		// TODO error occurs here
		velocity[i] -= dt * sum;
	}
}

float Fluid::calculateVmax() {
	float vmax = 0.f;
	for (int i = 0; i < numParticles; i++) {
		const float magnitude = velocity[i].magnitude();
		if (magnitude > vmax) vmax = magnitude;
	}
	return vmax;
}

float Fluid::calculateTimeStep() {
	const float vmaxLimit = 100.f;
	const float vmax = calculateVmax();

	const float vmaxSafe = std::min(vmaxLimit, 2.f * vmax);
	return 0.4 * smoothingLen / (vmaxSafe + 1E-6);
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
	// std::vector<int> neighbors = findNeighbors(particleIndex);
	
	for (auto j : neighbors[particleIndex]) {
		if (particleIndex == j) continue;
		const float Ai = field[particleIndex] / (density[particleIndex] * density[particleIndex]);
		const float Aj = field[j] / (density[j] * density[j]);
		Vector2f grad = smoothingGradient(position[particleIndex] - position[j], smoothingLen);
		sum += mass[j] * (Ai + Aj) * grad;
		if (std::isnan(sum.x) || std::isnan(sum.y)) {
			printf("Ai %f, Aj %f, grad (%f, %f)\n", Ai, Aj, grad.x, grad.y);
		}
	}
	return density[particleIndex] * sum;
}

float Fluid::divergence(int particleIndex, std::vector<Vector2f> field) {
	float sum = 0.f;
	// std::vector<int> neighbors = findNeighbors(particleIndex);

	for (auto j : neighbors[particleIndex]) {
		if (particleIndex == j) continue;
		const Vector2f Aij = field[particleIndex] - field[j];
		const Vector2f grad = smoothingGradient(position[particleIndex] - position[j], smoothingLen);
		sum += mass[j] * Aij.dot(grad);
	}
	return -1.f / density[particleIndex] * sum;
}

Vector2f Fluid::laplacian(int particleIndex, std::vector<Vector2f> field) {
	Vector2f sum = Vector2f(0.f, 0.f);
	// std::vector<int> neighbors = findNeighbors(particleIndex);

	for (auto j : neighbors[particleIndex]) {
		if (particleIndex == j) continue;
		const float volume = mass[j] / density[j];
		const Vector2f Aij = field[particleIndex] - field[j];
		const Vector2f Xij = position[particleIndex] - position[j];
		const Vector2f smoothGrad = smoothingGradient(Xij, smoothingLen);
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


float Fluid::smoothingKernel(float dist, float smoothingLength) {
	if (dist < 0 || dist > smoothingLength) return 0.f;
	// poly6 kernel
	const float coeff = 4 / (M_PI * pow(smoothingLength, 8.f));
	const float factor = smoothingLength * smoothingLength - dist * dist;
	return coeff * factor * factor * factor;
}

Vector2f Fluid::smoothingGradient(Vector2f r, float smoothingLength) {
	const float dist = r.magnitude();
	if (dist < 0 || dist > smoothingLength) return Vector2f(0.f, 0.f);
	// poly6 kernel gradient
	const float coeff = 4 / (M_PI * pow(smoothingLength, 8.f));
	const float factor = -6 * (smoothingLength * smoothingLength - dist * dist);
	return coeff * (factor * factor) * r;
}