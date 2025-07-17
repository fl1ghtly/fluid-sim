#include "ParticleSystem.h"

ParticleSystem::ParticleSystem(int count) : vertices(sf::PrimitiveType::Points, count)
{
}

 void ParticleSystem::update(std::vector<Vector2f> position, std::vector<Vector2f> velocity) {
	// const Vector2f maxVel = *std::max_element(velocity.begin(), velocity.end());
	// const float magnitude = maxVel.magnitude();
	const float magnitude = 20.f;
	for (int i = 0; i < position.size(); i++) {
		vertices[i].position = {position[i].x, position[i].y};
		const float t = std::clamp(velocity[i].magnitude() / magnitude, 0.f, 1.f);
		vertices[i].color = linearGradient(t);
	}
 }

 void ParticleSystem::draw(sf::RenderTarget& target, sf::RenderStates states) const {
	states.texture = nullptr;
	target.draw(vertices, states);
 }

 sf::Color ParticleSystem::linearGradient(float t) {
	const sf::Color start = sf::Color::Blue;
	const sf::Color mid = sf::Color::Green;
	const sf::Color end = sf::Color::Red;

	sf::Color sampled; 
	if (t <= 0.5) {
		sampled.r = (mid.r * t * 2.f) + start.r * (0.5 - t) * 2.f;
		sampled.g = (mid.g * t * 2.f) + start.g * (0.5 - t) * 2.f;
		sampled.b = (mid.b * t * 2.f) + start.b * (0.5 - t) * 2.f;
	} else {
		sampled.r = end.r * (t - 0.5f) * 2.f + mid.r * (1.f - t) * 2.f;
		sampled.g = end.g * (t - 0.5f) * 2.f + mid.g * (1.f - t) * 2.f;
		sampled.b = end.b * (t - 0.5f) * 2.f + mid.b * (1.f - t) * 2.f;
	}

	return sampled;
 }