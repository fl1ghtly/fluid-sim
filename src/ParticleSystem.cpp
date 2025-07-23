#include "ParticleSystem.h"

ParticleSystem::ParticleSystem(int count) : vertices(sf::PrimitiveType::Points, count)
{
}

 void ParticleSystem::update(
	std::vector<Vector2f> position, 
	std::vector<Vector2f> field, 
	float min, 
	float max, 
	std::vector<sf::Color> cmap) {
	for (int i = 0; i < position.size(); i++) {
		vertices[i].position = {position[i].x, position[i].y};
		vertices[i].color = getColorMap(field[i].magnitude(), min, max, cmap);
	}
 }

 void ParticleSystem::update(
	std::vector<Vector2f> position, 
	std::vector<float> field, 
	float min, 
	float max, 
	std::vector<sf::Color> cmap) {
	for (int i = 0; i < position.size(); i++) {
		vertices[i].position = {position[i].x, position[i].y};
		vertices[i].color = getColorMap(field[i], min, max, cmap);
	}
 }

 void ParticleSystem::draw(sf::RenderTarget& target, sf::RenderStates states) const {
	states.texture = nullptr;
	target.draw(vertices, states);
 }

 sf::Color ParticleSystem::getColorMap(float v, float vmin, float vmax, std::vector<sf::Color> cmap) {
	if (v <= vmin) return cmap[0];
	if (v >= vmax) return cmap.back();

	const float interval = (vmax - vmin) / (cmap.size() - 1);
	
	const int startIndex = (int) (v - vmin) / interval; 
	sf::Color c1 = cmap[startIndex];
	sf::Color c2 = cmap[startIndex + 1];
	const float t = (v - vmin - startIndex * interval) / interval;

	return sf::Color(
		(1.f - t) * c1.r + t * c2.r,
		(1.f - t) * c1.g + t * c2.g,
		(1.f - t) * c1.b + t * c2.b
	);
 }