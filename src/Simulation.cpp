#include "Simulation.h"

Simulation::Simulation(
	float width, 
	float height, 
	FluidParameters& params,
	Vector2f downDir, 
	float fixedTimestep
) 
	: width(width), 
	  height(height), 
	  currentStep(0),
	  numParticles(0),
	  params(params),
	  downDir(downDir),
	  fixedTimestep(fixedTimestep)
{
	constexpr float PI = std::numbers::pi_v<float>;
	// Kernel Coefficients
	poly6C = 4.f / (PI * pow(params.smoothingRadius, 8.f));
	spikyGC = -30.f / (PI * pow(params.smoothingRadius, 5.f));
	viscosityLC = 40.f / (PI * pow(params.smoothingRadius, 5.f));
	spatialGrid = createHashmap<GridCell, std::vector<std::size_t>>();
	boundarySpatialGrid = createHashmap<GridCell, std::vector<std::size_t>>();
}

void Simulation::initializeParticleValues(std::size_t amount) {
	const float particleArea = (4.f / 9.f) * params.smoothingRadius * params.smoothingRadius;
	const float particleMass = params.restDensity * particleArea;

	numParticles += amount;
	position.reserve(numParticles);
	velocity.resize(numParticles, {0.f, 0.f});
	mass.resize(numParticles, particleMass);
	density.resize(numParticles, params.restDensity);
	pressure.resize(numParticles, 0.f);
	neighbors.resize(numParticles);
	boundaryNeighbors.resize(numParticles);
}

void Simulation::initializeParticleGrid(Vector2f center, std::size_t gridWidth, std::size_t amount) {
	// Spacing factor of 2.5 results in about 20 neighbors on average
	constexpr float spacingFactor = 2.5f;
	const float spacing = params.smoothingRadius / spacingFactor;
	
	const float gridHeight = static_cast<float>(amount) / static_cast<float>(gridWidth);
	
	initializeParticleValues(amount);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> jitter(0, spacing);
	
	// X and Y alone define the top left corner of grid, therefore
	// we need to translate the center of the grid to the top left
	center.x -= spacing * static_cast<float>(gridWidth) / 2.f;
	center.y -= spacing * gridHeight / 2.f;
	for (std::size_t i = 0; i < amount; i++) {
		position.push_back(Vector2f(
			center.x + jitter(gen) + spacing * static_cast<float>(i % gridWidth),
			center.y + jitter(gen) + spacing * static_cast<float>(i / gridWidth)
		));
	}
}

void Simulation::initializeParticleRandom(Vector2f min, Vector2f max, std::size_t amount) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> xDistr(min.x, max.x);
	std::uniform_real_distribution<float> yDistr(min.y, max.y);

	initializeParticleValues(amount);
	
	for (std::size_t i = 0; i < amount; i++) {
		position.push_back(Vector2f(xDistr(gen), yDistr(gen)));
	}
}

void Simulation::initializeParticleCircle(Vector2f origin, float radius) {
	constexpr float spacingFactor = 1.5f;
	const float spacing = params.smoothingRadius / spacingFactor;
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> jitter(0, spacing);
	
	const float sqRadius = radius * radius;
	std::size_t amount = 0;
	// Calculate only for a quarter circle
	for (float x = 0; x < radius; x += spacing) {
		const float xSq = x * x;
		for (float y = 0; y < radius; y += spacing) {
			if (xSq + y * y > sqRadius) continue;
			const Vector2f point = {
				origin.x + jitter(gen),
				origin.y + jitter(gen)
			};
			position.push_back({point.x + x, point.y + y});
			position.push_back({point.x - x, point.y + y});
			position.push_back({point.x - x, point.y - y});
			position.push_back({point.x + x, point.y - y});
			amount += 4;
		}
	}
	initializeParticleValues(amount);
}

void Simulation::addBoundary(Boundary b) {
	const auto pos = b.getBoundaryParticlePositions();
	const auto vol = b.getBoundaryParticleVolume(); 

	boundaries.push_back(b);
	// Keep track of which indices belong to a boundary object by storing the
	// starting index for the boundary
	boundaryStartIndex.push_back(boundaryPos.size());
	boundaryPos.insert(boundaryPos.end(), pos.begin(), pos.end());
	boundaryVolume.insert(boundaryVolume.end(), vol.begin(), vol.end());
	boundaryForce.resize(boundaries.size());
}

void Simulation::addBoundary(std::vector<Boundary> b) {
	boundaries.insert(boundaries.end(), b.begin(), b.end());
	boundaryForce.resize(boundaries.size());
	for (const auto& boundary : b) {
		boundaryStartIndex.push_back(boundaryPos.size());
		const auto pos = boundary.getBoundaryParticlePositions();
		const auto vol = boundary.getBoundaryParticleVolume();
		boundaryPos.insert(boundaryPos.end(), pos.begin(), pos.end());
		boundaryVolume.insert(boundaryVolume.end(), vol.begin(), vol.end());
	}
}

void Simulation::update() {
	SimZoneScoped;
	if (numParticles <= 0) return;
	const float deltaTime = calculateTimeStep();

	findNeighbors();
	applyNonPressureForce(deltaTime);
	calculateDensity(deltaTime);
	calculatePressure();
	applyPressureForce(deltaTime);
	applyForceToBoundary(deltaTime);
	applyWorldBoundary();
	currentStep++;
}

void Simulation::calculateDensity(float dt) {
	SimZoneScoped;
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		density[i] = 0.f;
		
		// Density from fluid particles
		for (auto j : neighbors[i]) {
			const Vector2f rij = position[i] - position[j];
			const float dist = rij.magnitude();
			const float influence = poly6Kernel(dist, params.smoothingRadius);
			const Vector2f velDiff = velocity[i] - velocity[j];
			const Vector2f smoothGrad = poly6Gradient(rij, params.smoothingRadius);
			density[i] += mass[j] * influence + dt * velDiff.dot(smoothGrad);
		}

		// Density from boundary particles
		for (auto k : boundaryNeighbors[i]) {
			const Vector2f rik = position[i] - boundaryPos[k];
			const float dist = rik.magnitude();
			const float influence = poly6Kernel(dist, params.smoothingRadius);
			density[i] += params.restDensity * boundaryVolume[k] * influence;
		}
	}
}

void Simulation::calculatePressure() {
	SimZoneScoped;
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		pressure[i] = params.stiffness * (pow(density[i] / params.restDensity, 7.f) - 1);
	}
}

void Simulation::applyNonPressureForce(float dt) {
	SimZoneScoped;
	std::vector<Vector2f> forces(numParticles);
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		// Force due to gravity
		Vector2f gForce(0.f, 0.f); 
		gForce = params.gravity * downDir * mass[i];
		
		// Force due to viscosity
		Vector2f vForce(0.f, 0.f);
		for (auto j : neighbors[i]) {
			// TODO missing extra density
			const float volume = mass[j] / density[j];
			const Vector2f velDiff = velocity[j] - velocity[i];
			const float dist = (position[i] - position[j]).magnitude();
			const float vLaplacian = viscosityLaplacian(dist, params.smoothingRadius);
			vForce += volume * velDiff * vLaplacian;
		}
		vForce *= params.viscosity;

		// Force due to viscosity between fluid and boundary
		Vector2f vbForce(0.f, 0.f);
		for (auto k : boundaryNeighbors[i]) {
			const Vector2f velDiff = -1 * velocity[i];
			const float dist = (position[i] - boundaryPos[k]).magnitude();
			const float vLaplacian = viscosityLaplacian(dist, params.smoothingRadius);
			vbForce += boundaryVolume[k] * velDiff * vLaplacian;
		}
		vbForce *= params.viscosity;
		
		forces[i] = gForce + vForce + vbForce;
	}
	
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		velocity[i] += dt * forces[i] / mass[i];
		position[i] += dt * velocity[i];
	}
}

void Simulation::applyPressureForce(float dt) {
	SimZoneScoped;
	std::vector<Vector2f> forces(numParticles);

	// Force due to pressure from other fluids
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		Vector2f sum(0.f, 0.f);
		const float p_rho_i = pressure[i] / (density[i] * density[i]);
		for (auto j : neighbors[i]) {
			if (j == i) continue;
			const float p_rho_j = pressure[j] / (density[j] * density[j]);
			const float combined = p_rho_i + p_rho_j;
			const Vector2f rij = position[i] - position[j];
			const Vector2f grad = spikyGradient(rij, params.smoothingRadius);
			sum += mass[j] * combined * grad;
		}
		forces[i] = -mass[i] * sum;
	}

	// Force due to pressure from boundary
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		Vector2f sum(0.f, 0.f);
		const float p_rho = pressure[i] / (density[i] * density[i]);
		for (auto k : boundaryNeighbors[i]) {
			const Vector2f rik = position[i] - boundaryPos[k];
			const Vector2f grad = spikyGradient(rik, params.smoothingRadius);
			const float boundaryMass = params.restDensity * boundaryVolume[k];
			sum += boundaryMass * grad;
			// Get the iterator pointing to the first start index greater than k
			const auto boundaryUpper = std::upper_bound(
				boundaryStartIndex.begin(),
				boundaryStartIndex.end(),
				k
			);
			const auto prevIndex = static_cast<std::size_t>(
				std::distance(boundaryStartIndex.begin(), boundaryUpper) - 1
			);
			// Accumulate force for k's boundary object
			boundaryForce[prevIndex] += mass[i] * p_rho * boundaryMass * grad;
		}
		forces[i] -= mass[i] * p_rho * sum;
	}
	
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		velocity[i] += dt * forces[i] / mass[i];
		position[i] += dt * velocity[i];
	}
}

void Simulation::applyForceToBoundary(float dt) {
	for (std::size_t i = 0; i < boundaries.size(); i++) {
		boundaries[i].applyForceAndTorque(boundaryForce[i], dt);
	}
}

void Simulation::applyWorldBoundary() {
	SimZoneScoped;
	for (std::size_t i = 0; i < numParticles; i++) {
		auto& p = position[i];
		auto& v = velocity[i];

		if (p.x <= 0.f) {
			p.x = 0.f;
			v.x *= -params.damping;
		} else if (p.x >= width) {
			p.x = width - 1.f;
			v.x *= -params.damping;
		}
		
		if (p.y <= 0.f) {
			p.y = 0.f;
			v.y *= -params.damping;
		} else if (p.y >= height) {
			p.y = height - 1.f;
			v.y *= -params.damping;
		}
	}
}

float Simulation::calculateTimeStep() {
	if (fixedTimestep > 0.f) return fixedTimestep;
	const Vector2f v = *std::max_element(velocity.begin(), velocity.end());
	const float vMax = v.magnitude();
	const float timestep = 0.4f * 2.f * params.smoothingRadius / (vMax + 1E-6f);

	return std::clamp(timestep, 1E-6f, 1.f/120.f);
}

void Simulation::buildSpatialGrid() {
	SimZoneScoped;
	// Cell Size = Kernel Support Length is optimal
	const float cellSize = params.smoothingRadius;

	// Build spatial grid for fluid particles
	spatialGrid->clear();
	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		GridCell c = {position[i], cellSize};
		if (spatialGrid->find(c)) {
			#pragma omp critical
			spatialGrid->get(c).push_back(i);
		} else {
			spatialGrid->insert(c, {i});
		}
	}
	
	// Build spatial grid for boundary particles
	boundarySpatialGrid->clear();
	#pragma omp parallel for
	for (std::size_t i = 0; i < boundaryPos.size(); i++) {
		GridCell c = {boundaryPos[i], cellSize};
		if (boundarySpatialGrid->find(c)) {
			#pragma omp critical
			boundarySpatialGrid->get(c).push_back(i);
		} else {
			boundarySpatialGrid->insert(c, {i});
		}
	}
}

void Simulation::findNeighbors() {
	SimZoneScoped;
	const float cellSize = params.smoothingRadius;
	// Maintain spatial locality every simulation step
	sortZIndex();
	buildSpatialGrid();

	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
		std::vector<std::size_t> neighborIndices;
		std::vector<std::size_t> boundaryNeighborIndices;
		// Reserve a bit more than the expected average number of neighbors (20)
		// to avoid constant reallocations
		neighborIndices.reserve(64);
		boundaryNeighborIndices.reserve(64);
		const GridCell center = {position[i], cellSize};
		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				const GridCell cell = {
					static_cast<int>(center.x) + dx, 
					static_cast<int>(center.y) + dy, 
					cellSize
				};
				const auto cellIndices = spatialGrid->get(cell);
				const auto boundaryCellIndices = boundarySpatialGrid->get(cell);
				
				for (const auto neighbor : cellIndices) {
					const Vector2f rij = position[i] - position[neighbor];
					const float sqDist = rij.dot(rij);
					if (sqDist > params.sqSmoothingRadius) continue;
					neighborIndices.push_back(neighbor);
				}

				for (const auto neighbor : boundaryCellIndices) {
					const Vector2f rij = position[i] - boundaryPos[neighbor];
					const float sqDist = rij.dot(rij);
					if (sqDist > params.sqSmoothingRadius) continue;
					boundaryNeighborIndices.push_back(neighbor);
				}
			}
		}
		neighbors[i] = std::move(neighborIndices);
		boundaryNeighbors[i] = std::move(boundaryNeighborIndices);
	}
}

void Simulation::sortZIndex() {
	SimZoneScoped;
	const float cellSize = params.smoothingRadius;
	std::vector<std::size_t> zIndex(numParticles);

	#pragma omp parallel for
	for (std::size_t i = 0; i < numParticles; i++) {
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

std::vector<Vector2f> Simulation::getPosition() {
	return position;
}

std::vector<Vector2f> Simulation::getVelocity() {
	return velocity;
}

std::vector<float> Simulation::getPressure() {
	return pressure;
}

std::vector<float> Simulation::getDensity() {
	return density;
}

std::size_t Simulation::getNumParticles() {
	return numParticles;
}

void Simulation::setDownDirection(Vector2f d) {
	downDir = d.normalize();
}

float Simulation::getPressureAtPoint(const Vector2f point) {
	const float cellSize = params.smoothingRadius;
	const GridCell center = {point, cellSize};
	float p = 0.f;

	for (int dx = -1; dx <= 1; dx++) {
		for (int dy = -1; dy <= 1; dy++) {
			GridCell cell = {
				static_cast<int>(center.x) + dx, 
				static_cast<int>(center.y) + dy, 
				cellSize
			};
			auto cellIndices = spatialGrid->get(cell);

			for (auto neighbor : cellIndices) {
				p += pressure[neighbor];
			}
		}
	}
	return p;
}

float Simulation::poly6Kernel(float dist, float smoothingLength) {
	if (dist < 0 || dist > smoothingLength) return 0.f;
	const float factor = params.sqSmoothingRadius - dist * dist;
	return poly6C * factor * factor * factor;
}

Vector2f Simulation::poly6Gradient(Vector2f r, float smoothingLength) {
	const float sqDist = r.dot(r);
	const float sqLength = smoothingLength * smoothingLength;
	if (sqDist < 0 || sqDist > sqLength) return Vector2f(0.f, 0.f);
	const float factor = -6 * (sqLength - sqDist);
	return poly6C * factor * factor * r;
}

Vector2f Simulation::spikyGradient(Vector2f r, float smoothingLength) {
	const float sqDist = r.dot(r);
	if (sqDist <= 0 || sqDist > params.sqSmoothingRadius) return Vector2f(0.f, 0.f);
	const float dist = sqrt(sqDist);
	const float factor = smoothingLength - dist;
	return spikyGC * factor * factor * r / dist;
}

float Simulation::viscosityLaplacian(float dist, float smoothingLength) {
	if (dist < 0 || dist > smoothingLength) return 0.f;
	const float factor = smoothingLength - dist;
	return viscosityLC * factor;
}