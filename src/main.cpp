#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "fluid.h"
#include "Vector2f.h"
#include "ParticleSystem.h"

int main(void) {
    constexpr int width = 200;
    constexpr int height = 600;
    constexpr float radius = 8.f;
    constexpr float fluidDensity = 1000.f;
    constexpr int numParticles = 1000;

	auto window = sf::RenderWindow(sf::VideoMode({width, height}), "Fluid Simulation");
    window.setFramerateLimit(144);
    
    ParticleSystem particles(numParticles);
	Fluid fluid(width, height, numParticles, radius, fluidDensity);

    // Create square grids
	fluid.initializeParticleGrid((int)sqrt(numParticles));
    // fluid.initializeParticleRandom();

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
		std::vector<Vector2f> vel = fluid.getVelocity();
        particles.update(pos, vel);
        
        window.clear();
        window.draw(particles);
        window.display();
        fluid.update();
    }
}