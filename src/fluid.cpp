#include "fluid.h"

#define GRAVITY 9.8
#define DAMPING 0.5

constexpr Vector2f VecDown = Vector2f(0.f, 1.f);
constexpr float h = 1.5f;

Fluid::Fluid(int width, int height, int numParticles, float radius) 
	: width(width), height(height), numParticles(numParticles),
	  radius(radius), position(numParticles), velocity(numParticles),
	  mass(numParticles), density(numParticles), pressure(numParticles)
{
	oldTime = clock();
}

void Fluid::initializeParticleValues() {
	std::fill(velocity.begin(), velocity.end(), Vector2f(0.f, 0.f));
	std::fill(mass.begin(), mass.end(), 1000.f);
	std::fill(density.begin(), density.end(), 0.f);
	std::fill(pressure.begin(), pressure.end(), 0.f);
	calculateDensity();
	calculatePressure();
}

void Fluid::initializeParticleGrid(int x, int y, int spacing, int width, int height) {
	for (int i = 0; i < numParticles; i++) {
		position[i] = Vector2f(x + spacing * (i % width), y + spacing * (i / height));
	}
	initializeParticleValues();
}

void Fluid::update() {
	clock_t currentTime = clock();
	float deltaTime = (float)(currentTime - oldTime) / CLOCKS_PER_SEC;
	std::vector<float> alpha = calculateAlphaFactors();
	applyNonPressureForce(deltaTime);
	// TODO adapt timestep to CFL condition
	correctDensityError(alpha, deltaTime);
	// Update Positions
	for (int i = 0; i < numParticles; i++) {
		position[i] += deltaTime * velocity[i];
	}
	calculateDensity();
	alpha = calculateAlphaFactors();
	correctDivergenceError(alpha, deltaTime);

	for (int i = 0; i < numParticles; i++) {
		// auto& v = velocity[i];
		auto& p = position[i];
		
		// v += Vector2f(0.f, 1.f) * GRAVITY * deltaTime;
		// p += v * deltaTime;

		if (p.x < 0) {
			p.x = 0;
		} else if (p.x >= width) {
			p.x = width - 1;
			// v.x *= -1 * DAMPING;
		}

		if (p.y < 0) {
			p.y = 0;
		} else if (p.y >= height) {
			p.y = height - 1;
			// v.y *= -1 * DAMPING;
		}
	}
	oldTime = currentTime;
}

std::vector<float> Fluid::calculateAlphaFactors() {
	std::vector<float> alpha(numParticles);
	for (int i = 0; i < numParticles; i++) {
		Vector2f sumLeft = Vector2f(0.f, 0.f);
		float sumRight = 0.f;
		for (int j = 0; j < numParticles; j++) {
			if (i == j) continue;
			const Vector2f mWij = mass[j] * smoothingGradient(position[i] - position[j], h);
			sumLeft += mWij;
			sumRight += mWij.magnitude() * mWij.magnitude();
		}
		const float lhs = sumLeft.magnitude() * sumLeft.magnitude();
		alpha[i] = std::max((float)10E-6, density[i] / (lhs + sumRight));
	}
	return alpha;
}

void Fluid::calculateDensity() {
	// TODO optimize
	for (int i = 0; i < numParticles; i++) {
		density[i] = 0.f;
		for (int j = 0; j < numParticles; j++) {
			if (i == j) continue;
			float dist = (position[i] - position[j]).magnitude();
			float influence = smoothingKernel(dist, h);
			density[i] += mass[j] * influence;
		}
	}
}

void Fluid::calculatePressure() {
	constexpr float restDensity = 1000.f;
	constexpr float stiffness = 10.f;
	for (int i = 0; i < numParticles; i++) {
		// Rest density for water is 1000 kg/m^3
		pressure[i] = stiffness * (density[i] - restDensity);
	}
}

void Fluid::applyNonPressureForce(float dt) {
	for (int i = 0; i < numParticles; i++) {
		// Force due to gravity
		Vector2f Force = VecDown * GRAVITY * mass[i];
		// TODO add other forces (Viscosity)
	
		velocity[i] += dt * Force / mass[i];
	}
}

void Fluid::correctDivergenceError(std::vector<float> alpha, float dt) {
	std::vector<float> densityChange;
	// TODO loop until divergence error less than threshold
	for (int i = 0; i < numParticles; i++) {
		densityChange.push_back(-density[i] * divergence(i, velocity));
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
			sum += mass[j] * stiffIJ * smoothingGradient(position[i] - position[j], h);
		}
		velocity[i] -= dt * sum;
	}
}

void Fluid::correctDensityError(std::vector<float> alpha, float dt) {
	// TODO loop until density error greater than some threshold
	const float restDensity = 1000.f;
	calculateDensity();
	for (int i = 0; i < numParticles; i++) {
		Vector2f sum = Vector2f(0.f, 0.f);
		const float stiffnessI = alpha[i] * (density[i] - restDensity) / (dt * dt);
		for (int j = 0; j < numParticles; j++) {
			if (i == j) continue;
			const float stiffnessJ = alpha[j] * (density[j] - restDensity) / (dt * dt);
			const float stiffIJ = stiffnessI / density[i] + stiffnessJ / density[j];
			sum += mass[j] * stiffIJ * smoothingGradient(position[i] - position[j], h);
		}
		// TODO error occurs here
		velocity[i] -= dt * sum;
	}
}

template <typename T>
Vector2f Fluid::gradient(int particleIndex, std::vector<T> field) {
	Vector2f sum = Vector2f(0.f, 0.f);
	for (int j = 0; j < numParticles; j++) {
		if (particleIndex == j) continue;
		const float Ai = field[particleIndex] / (density[particleIndex] * density[particleIndex]);
		const float Aj = field[j] / (density[j] * density[j]);
		sum += mass[j] * (Ai + Aj) * smoothingGradient(position[particleIndex] - position[j], h);
	}
	return density[particleIndex] * sum;
}

float Fluid::divergence(int particleIndex, std::vector<Vector2f> field) {
	float sum = 0.f;
	for (int j = 0; j < numParticles; j++) {
		if (particleIndex == j) continue;
		const Vector2f Aij = field[particleIndex] - field[j];
		const Vector2f grad = smoothingGradient(position[particleIndex] - position[j], h);
		sum += mass[j] * Aij.dot(grad);
	}
	return -1.f / density[particleIndex] * sum;
}

std::vector<Vector2f> Fluid::getPosition() {
	return position;
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