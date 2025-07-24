#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <algorithm>
#include "Vector2f.h"

class ParticleSystem : public sf::Drawable {
	public:
		ParticleSystem();
		void update(std::vector<Vector2f> position, std::vector<Vector2f> field, float min, float max, std::vector<sf::Color> cmap);
		void update(std::vector<Vector2f> position, std::vector<float> field, float min, float max, std::vector<sf::Color> cmap);
		void update(std::vector<Vector2f> position, sf::Color color);
		
	private:
		sf::VertexArray vertices;

		void draw(sf::RenderTarget& target, sf::RenderStates states) const override;
		sf::Color getColorMap(float v, float vmin, float vmax, std::vector<sf::Color> cmap);
};