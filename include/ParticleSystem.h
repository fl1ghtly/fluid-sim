#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <algorithm>
#include "Vector2f.h"

class ParticleSystem : public sf::Drawable {
	public:
		ParticleSystem(int count);
		void update(std::vector<Vector2f> position, std::vector<Vector2f> velocity);
		
	private:
		sf::VertexArray vertices;

		void draw(sf::RenderTarget& target, sf::RenderStates states) const override;
		sf::Color linearGradient(float t);
};