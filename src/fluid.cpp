#include "fluid.h"

Fluid::Fluid(int width, int height, FluidParameters& params, float fixedTimestep) 
	: width(width), 
	  height(height), 
	  params(params),
	  fixedTimestep(fixedTimestep),
	  position(params.numParticles), 
	  velocity(params.numParticles), 
	  mass(params.numParticles), 
	  density(params.numParticles), 
	  pressure(params.numParticles),
	  neighbors(params.numParticles)
{
	// Kernel Coefficients
	poly6C = 4 / (M_PI * pow(params.smoothingRadius, 8.f));
	spikyGC = -30.f / (M_PI * pow(params.smoothingRadius, 5.f));
	viscosityLC = 40.f / (M_PI * pow(params.smoothingRadius, 5.f));
}

void Fluid::initializeParticleValues() {
	const float particleArea = (4.f / 9.f) * params.smoothingRadius * params.smoothingRadius;
	const float particleMass = params.restDensity * particleArea;
	std::fill(velocity.begin(), velocity.end(), Vector2f(0.f, 0.f));
	std::fill(mass.begin(), mass.end(), particleMass);
	std::fill(density.begin(), density.end(), params.restDensity);
	std::fill(pressure.begin(), pressure.end(), 0.f);
}

void Fluid::initializeParticleGrid(int gridWidth) {
	// Spacing factor of 2.5 results in about 20 neighbors on average
	constexpr float spacingFactor = 2.5f;
	const float spacing = params.smoothingRadius / spacingFactor;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> jitter(0, spacing);
	
	Vector2f center(width / 2, height / 2);
	center.x -= 2 * gridWidth;
	center.y -= 2 * params.numParticles / gridWidth;
	for (int i = 0; i < params.numParticles; i++) {
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
	for (int i = 0; i < params.numParticles; i++) {
		position[i] = Vector2f(xDistr(gen), yDistr(gen));
	}
	initializeParticleValues();
}

void Fluid::update() {
	const float deltaTime = calculateTimeStep();
	findNeighbors();
	applyNonPressureForce(deltaTime);
	calculateDensity(deltaTime);
	calculatePressure();
	applyPressureForce(deltaTime);       
	applyBoundaryCondition();
}

void Fluid::calculateDensity(float dt) {
	ZoneScoped;
	#pragma omp parallel for
	for (int i = 0; i < params.numParticles; i++) {
		density[i] = 0.f;
		for (auto j : neighbors[i]) {
			const Vector2f rij = position[i] - position[j];
			const float dist = rij.magnitude();
			const float influence = poly6Kernel(dist, params.smoothingRadius);
			const Vector2f velDiff = velocity[i] - velocity[j];
			const Vector2f smoothGrad = poly6Gradient(rij, params.smoothingRadius);
			density[i] += mass[j] * influence + dt * velDiff.dot(smoothGrad);
		}
	}
}

void Fluid::calculatePressure() {
	ZoneScoped;
	#pragma omp parallel for
	for (int i = 0; i < params.numParticles; i++) {
		pressure[i] = params.stiffness * (pow(density[i] / params.restDensity, 7.f) - 1);
	}
}

void Fluid::applyNonPressureForce(float dt) {
	ZoneScoped;
	std::vector<Vector2f> forces(params.numParticles);
	#pragma omp parallel for
	for (int i = 0; i < params.numParticles; i++) {
		// Force due to gravity
		Vector2f gForce(0.f, 0.f); 
		gForce = params.getGravityVector() * mass[i];
		
		// Force due to viscosity
		Vector2f vForce(0.f, 0.f);
		for (int j : neighbors[i]) {
			const float volume = mass[j] / density[j];
			const Vector2f velDiff = velocity[j] - velocity[i];
			const float dist = (position[i] - position[j]).magnitude();
			const float vLaplacian = viscosityLaplacian(dist, params.smoothingRadius);
			vForce += volume * velDiff * vLaplacian;
		}
		vForce *= params.viscosity;
		
		forces[i] = gForce + vForce;
	}
	
	#pragma omp parallel for
	for (int i = 0; i < params.numParticles; i++) {
		velocity[i] += dt * forces[i] / mass[i];
		position[i] += dt * velocity[i];
	}
}

void Fluid::applyPressureForce(float dt) {
	std::vector<Vector2f> forces(params.numParticles);
	#pragma omp parallel for
	for (int i = 0; i < params.numParticles; i++) {
		Vector2f force(0.f, 0.f);
		for (int j : neighbors[i]) {
			const float volume = mass[j] / density[j];
			const float avgPressure = (pressure[i] + pressure[j]) / 2.f;
			const Vector2f rij = position[i] - position[j];
			const Vector2f grad = spikyGradient(rij, params.smoothingRadius);
			force += -volume * avgPressure * grad;
		}
		forces[i] =  force;
	}
	
	#pragma omp parallel for
	for (int i = 0; i < params.numParticles; i++) {
		velocity[i] += dt * forces[i] / mass[i];
		position[i] += dt * velocity[i];
	}
}

void Fluid::applyBoundaryCondition() {
	for (int i = 0; i < params.numParticles; i++) {
		auto& p = position[i];
		auto& v = velocity[i];

		if (p.x <= 0.f) {
			p.x = 0.f;
			v.x *= -params.damping;
		} else if (p.x >= width) {
			p.x = width - 1;
			v.x *= -params.damping;
		}
		
		if (p.y <= 0.f) {
			p.y = 0.f;
			v.y *= -params.damping;
		} else if (p.y >= height) {
			p.y = height - 1;
			v.y *= -params.damping;
		}
	}
}

float Fluid::calculateTimeStep() {
	if (fixedTimestep > 0.f) return fixedTimestep;
	const Vector2f v = *std::max_element(velocity.begin(), velocity.end());
	const float vMax = v.magnitude();
	const float timestep = 0.4f * 2.f * params.smoothingRadius / (vMax + 1E-6f);

	return std::clamp(timestep, 1E-6f, 1.f/120.f);
}

void Fluid::buildSpatialGrid() {
	// Cell Size = Kernel Support Length is optimal
	const float cellSize = params.smoothingRadius;

	spatialGrid.clear();
	#pragma omp parallel for
	for (int i = 0; i < params.numParticles; i++) {
		GridCell c = {position[i], cellSize};
		spatialGrid[c].push_back(i);
	}
}

void Fluid::findNeighbors() {
	const float cellSize = params.smoothingRadius;
	buildSpatialGrid();
	
	for (int i = 0; i < params.numParticles; i++) {
		std::vector<int> neighborIndices;
		const GridCell center = {position[i], cellSize};
		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				GridCell cell = {center.x + dx, center.y + dy, cellSize};
				const std::vector<int> cellIndices = spatialGrid[cell];
				for (auto neighbor : cellIndices) {
					// if (neighbor == i) continue;
					const float dist = (position[i] - position[neighbor]).magnitude();
					if (dist > params.smoothingRadius) continue;
					neighborIndices.push_back(neighbor);
				}
			}
		}
		neighbors[i] = neighborIndices;
	}
	}
}

void Fluid::sortZIndex() {
	ZoneScoped;
	const float cellSize = params.smoothingRadius;
	std::vector<unsigned int> zIndex(params.numParticles);

	for (int i = 0; i < params.numParticles; i++) {
		const GridCell c = {position[i], cellSize};
		zIndex[i] = c.zOrder;
	}

	auto zipped = std::views::zip(
		zIndex, 
		position, 
		velocity,
		mass,
		density,
		pressure 
	);
	
	std::ranges::sort(zipped, [](const auto& lhs, const auto& rhs) {
		// Sort particle arrays based on their Z-index
		return std::get<0>(lhs) < std::get<0>(rhs);
	});
}

template <typename T>
Vector2f Fluid::gradient(int particleIndex, std::vector<T> field) {
	Vector2f sum = Vector2f(0.f, 0.f);
	
	for (auto j : neighbors[particleIndex]) {
		if (particleIndex == j) continue;
		const float Ai = field[particleIndex] / (density[particleIndex] * density[particleIndex]);
		const float Aj = field[j] / (density[j] * density[j]);
		Vector2f grad = spikyGradient(position[particleIndex] - position[j], params.smoothingRadius);
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
		const Vector2f grad = poly6Gradient(position[particleIndex] - position[j], params.smoothingRadius);
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
		const Vector2f smoothGrad = poly6Gradient(Xij, params.smoothingRadius);
		const float numerator = Xij.dot(smoothGrad);
		const float denominator = Xij.dot(Xij) + 0.01 * params.smoothingRadius * params.smoothingRadius;
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
	const float cellSize = params.smoothingRadius;
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
	const float factor = params.sqSmoothingRadius - dist * dist;
	return poly6C * factor * factor * factor;
}

Vector2f Fluid::poly6Gradient(Vector2f r, float smoothingLength) {
	const float sqDist = r.dot(r);
	if (sqDist <= 0 || sqDist > params.sqSmoothingRadius) return Vector2f(0.f, 0.f);
	const float factor = -6 * (params.sqSmoothingRadius - sqDist);
	return poly6C * factor * factor * r;
}

Vector2f Fluid::spikyGradient(Vector2f r, float smoothingLength) {
	const float sqDist = r.dot(r);
	if (sqDist <= 0 || sqDist > params.sqSmoothingRadius) return Vector2f(0.f, 0.f);
	const float dist = sqrt(sqDist);
	const float factor = smoothingLength - dist;
	return spikyGC * factor * factor * r / dist;
}

float Fluid::viscosityLaplacian(float dist, float smoothingLength) {
	if (dist <= 0 || dist > smoothingLength) return 0.f;
	const float factor = smoothingLength - dist;
	return viscosityLC * factor;
}