#include "PressureSystem.h"

PressureSystem::PressureSystem(int w, int h, int gridSize) 
	: vertices(sf::PrimitiveType::TriangleStrip),
	  width(w),
	  height(h),
	  gridSize(gridSize) 
{
	gridX = (float) width / gridSize;
	gridY = (float) height / gridSize;
}

void PressureSystem::draw(sf::RenderTarget& target, sf::RenderStates states) const {
	states.texture = nullptr;
	target.draw(vertices, states);
}

void PressureSystem::updatePressureGradient(std::vector<std::vector<float>> pressureGrid, float max, float min) {
	sf::Vector2f lastVertex;
	for (int y = 0; y < gridSize; y++) {
        // Add degenerate triangle (except first row)
        if (y > 0) {
            // Repeat last vertex of previous row and first vertex of new row
            vertices.append(sf::Vertex(lastVertex, sf::Color::White));
            vertices.append(sf::Vertex(sf::Vector2f(0, y * gridY), sf::Color::White));
        }

        // Add vertices for current and next row
        for (int x = 0; x <= gridSize; x++) {
            vertices.append(sf::Vertex(
                sf::Vector2f(x * gridX, y * gridY), 
                pressureToColor(pressureGrid[x][y], min, max)
            ));
			lastVertex = sf::Vector2f(x * gridX, (y + 1) * gridY);
            vertices.append(sf::Vertex(
				lastVertex, 
				pressureToColor(pressureGrid[x][y], min, max)
            ));
        }
    }
}

sf::Color PressureSystem::pressureToColor(const float p, const float min, const float max) {
	constexpr sf::Color start(6, 0, 85);
	constexpr sf::Color mid(255, 255, 255);
	constexpr sf::Color end(165, 7, 7);
	
	const float t = std::clamp((p - min) / (max - min), 0.f, 1.f);

	sf::Color sampled = mid; 
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