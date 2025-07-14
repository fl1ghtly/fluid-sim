#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "fluid.h"
#include "Vector2f.h"

int main(void) {
	auto window = sf::RenderWindow(sf::VideoMode({128u, 64u}), "CMake SFML Project");
    window.setFramerateLimit(144);
	
    constexpr float radius = 5.f;
    constexpr float fluidDensity = 1000.f;
    constexpr int numParticles = 100;

	Fluid fluid(128, 64, numParticles, radius, fluidDensity);
	fluid.initializeParticleGrid(0, 0, 20, 5);
    while (window.isOpen())
    {
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        
		std::vector<Vector2f> pos = fluid.getPosition();
        
        window.clear();
        
        for (auto p : pos) {
            sf::CircleShape shape(radius);
            shape.setFillColor(sf::Color(0, 0, 255));
            shape.setPosition({p.x, p.y});
            window.draw(shape);
        }
        
        window.display();
        fluid.update();
    }
}