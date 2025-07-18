#pragma once
#include <SFML/Graphics.hpp>
#include <vector>
#include <algorithm>

class PressureSystem : public sf::Drawable {
	public:
		PressureSystem(int width, int height, int gridSize);
		void updatePressureGradient(std::vector<std::vector<float>> pressureGrid, float max, float min);
		
	private:
		sf::VertexArray vertices;
		int width, height, gridSize;
		float gridX, gridY;

		void draw(sf::RenderTarget& target, sf::RenderStates states) const override;
		sf::Color pressureToColor(const float p, const float min, const float max);
};