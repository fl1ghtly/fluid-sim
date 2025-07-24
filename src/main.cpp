#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "Simulation.h"
#include "FluidParameters.h"
#include "Vector2f.h"
#include "ParticleSystem.h"
#include "PressureSystem.h"
#include "ColorMap.h"
#include "Boundary.h"

#include <tracy/Tracy.hpp>

int main(void) {
    constexpr Vector2f downDir(0.f, 1.f);
    constexpr float g = 9.8f;
    constexpr float damping = 0.95f;
    constexpr float restDensity = 1000.f;
    constexpr float stiffness = 1E+3f;
    constexpr float viscosity = 1E+6f;
    constexpr float smoothingRadius = 8.f;
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

    /*
    constexpr int GRID_SIZE = 50;
    constexpr float GRID_X = width / GRID_SIZE;
    constexpr float GRID_Y = height / GRID_SIZE;
    std::vector<std::vector<float>> pressureGrid;
    pressureGrid.resize(GRID_SIZE + 1, std::vector<float>(GRID_SIZE + 1, 0.f));
    PressureSystem pressureGradient(width, height, GRID_SIZE);
    */

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
    Boundary b1(10000, smoothingRadius);
    b1.createBox({0.25f * width, 0.6f * height}, {0.33f * width, 0.8f * height});

    // Triangle
    Boundary b2(10000, smoothingRadius);
    b2.createPolygon({{200.f, 800.f}, {400.f, 800.f}, {300.f, 1000.f}});

    // Circle
    Boundary b3(10000, smoothingRadius);
    b3.createCircle({width / 2.f, height - 200.f}, 100.f);

    std::vector<Boundary> boundaries = {b1, b2, b3};
    
    std::vector<Vector2f> boundaryPositions;
    for (const auto b : boundaries) {
        const auto pos = b.getBoundaryParticlePositions();
        boundaryPositions.insert(boundaryPositions.end(), pos.begin(), pos.end());
    }
    boundaryParticles.update(boundaryPositions, sf::Color::White);

    sim.addBoundary(boundaries);

    constexpr int MAX_SIM_STEPS = 1200;
    int steps = 0;

    while (window.isOpen())
    {
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        
		std::vector<Vector2f> pos = sim.getPosition();
		std::vector<Vector2f> vel = sim.getVelocity();
        particles.update(pos, vel, 0.f, 100.f, ColorMap::viridis);

        /*
        float maxP = 0.f;
        float minP = 0.f;
        for (int i = 0; i < GRID_SIZE + 1; i++) {
            for (int j = 0; j < GRID_SIZE + 1; j++) {
                const float p = sim.getPressureAtPoint({i * GRID_X, j * GRID_Y});
                pressureGrid[i][j] = p;
                if (p > maxP) maxP = p;
                if (p < minP) minP = p;
            }
        }
        pressureGradient.updatePressureGradient(pressureGrid, maxP, minP);
        */
        
        
        window.clear();
        // window.draw(pressureGradient);
        window.draw(particles);
        window.draw(boundaryParticles);
        window.display();
        sim.update();
        FrameMark;

        if (steps >= MAX_SIM_STEPS) {
            // window.close();
        }
        steps++;
    }
}