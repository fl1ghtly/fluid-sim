#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "fluid.h"
#include "FluidParameters.h"
#include "Vector2f.h"
#include "ParticleSystem.h"

int main(void) {
    constexpr Vector2f downDir(0.f, 1.f);
    constexpr float g = 9.8f;
    constexpr float damping = 0.95f;
    constexpr float restDensity = 1000.f;
    constexpr float stiffness = 10000.f;
    constexpr float viscosity = 1E+6f;
    constexpr float smoothingRadius = 8.f;
    constexpr int numParticles = 1000;
    constexpr int width = 800;
    constexpr int height = 600;

    FluidParameters params(
        downDir, 
        g, 
        damping, 
        restDensity, 
        stiffness, 
        viscosity, 
        smoothingRadius, 
        numParticles
    );

	auto window = sf::RenderWindow(sf::VideoMode({width, height}), "Fluid Simulation");
    window.setFramerateLimit(144);
    
    ParticleSystem particles(numParticles);
	Fluid fluid(width, height, params);

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