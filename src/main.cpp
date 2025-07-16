#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "fluid.h"
#include "Vector2f.h"
#include "ParticleSystem.h"

int main(void) {
    constexpr int width = 1920;
    constexpr int height = 1080;
    constexpr float radius = 1.2f;
    constexpr float fluidDensity = 1000.f;
    constexpr int numParticles = 2000;

	auto window = sf::RenderWindow(sf::VideoMode({width, height}), "Fluid Simulation");
    window.setFramerateLimit(144);
    
    ParticleSystem particles(numParticles);
	Fluid fluid(width, height, numParticles, radius, fluidDensity);
	fluid.initializeParticleGrid(width / 4, height / 4, 200, 5);
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