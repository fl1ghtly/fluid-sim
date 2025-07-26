#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "Simulation.h"
#include "FluidParameters.h"
#include "Vector2f.h"
#include "Visualization/ParticleSystem.h"
#include "Visualization/ColorMap.h"
#include "Boundary.h"

#ifdef TRACY_ENABLED
#include <tracy/Tracy.hpp>
#define SimFrameMark FrameMark

#else
#define SimFrameMark

#endif

int main(void) {
    constexpr Vector2f downDir(0.f, 1.f);
    constexpr float g = 9.8f;
    constexpr float damping = 0.95f;
    constexpr float restDensity = 1000.f;
    constexpr float stiffness = 1E+3f;
    constexpr float viscosity = 1E+6f;
    constexpr float smoothingRadius = 8.f;
    constexpr float boundaryRadius = smoothingRadius * 4.f;
    constexpr int numParticles = 50000;
    constexpr int width = 1920;
    constexpr int height = 1080;

    FluidParameters params(
        g, 
        damping, 
        restDensity, 
        stiffness, 
        viscosity, 
        smoothingRadius
    );

	auto window = sf::RenderWindow(sf::VideoMode({width, height}), "Simulation Simulation");
    window.setFramerateLimit(144);
    
	Simulation sim(width, height, params, downDir, -1.f);
    ParticleSystem particles;
    ParticleSystem boundaryParticles;
    
    // Create square grids
    sim.initializeParticleGrid(
        {width * 0.25f, height * 0.25f}, 
        (int)sqrt(numParticles / 2.f), 
        numParticles / 2.f
    );
    sim.initializeParticleGrid(
        {width * 0.75f, height * 0.75f}, 
        (int)sqrt(numParticles / 2.f), 
        numParticles / 2.f
    );
    
    // Box
    Boundary b1(boundaryRadius);
    b1.createBox({0.25f * width, 0.6f * height}, {0.33f * width, 0.8f * height}, 0.25);

    // Triangle
    Boundary b2(boundaryRadius);
    b2.createPolygon({{200.f, 800.f}, {400.f, 800.f}, {300.f, 1000.f}}, 0.0625);
    
    // Circle
    Boundary b3(boundaryRadius);
    b3.createCircle({width / 2.f, height - 200.f}, 100.f);
    
    std::vector<Boundary> boundaries = {b1, b2, b3};
    
    std::vector<Vector2f> boundaryPositions;
    for (auto& b : boundaries) {
        b.activateBoundary();
        const auto pos = b.getBoundaryParticlePositions();
        boundaryPositions.insert(boundaryPositions.end(), pos.begin(), pos.end());
    }
    boundaryParticles.update(boundaryPositions, sf::Color::White);
    
    sim.addBoundary(boundaries);

    while (window.isOpen()) {
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) window.close();
        }
        
        // Render fluid particles
		std::vector<Vector2f> pos = sim.getPosition();
		std::vector<Vector2f> vel = sim.getVelocity();
        particles.update(pos, vel, 0.f, 100.f, ColorMap::viridis);
        
        window.clear();
        window.draw(particles);
        window.draw(boundaryParticles);
        window.display();
        sim.update();
        SimFrameMark;
    }
}