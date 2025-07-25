#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "Simulation.h"
#include "FluidParameters.h"
#include "Vector2f.h"
#include "ParticleSystem.h"
#include "PressureSystem.h"
#include "ColorMap.h"
#include "Boundary.h"

#ifdef TRACY_ENABLED
#include <tracy/Tracy.hpp>
#define SimFrameMark FrameMark

#else
#define SimFrameMark

#endif

enum State {
    BOUNDARY,
    FLUID
};

enum BoundaryState {
    POLYGON,
    BOX,
    CIRCLE
};

struct BoundaryData {
    std::vector<Vector2f>& boundaryVertices;
    float smoothingRadius;
    float compression;
};

void handleFluidState(std::optional<sf::Event> event, Simulation& sim);
void handleBoundaryState(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
);

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
    
    /*
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
    for (const auto b : boundaries) {
        const auto pos = b.getBoundaryParticlePositions();
        boundaryPositions.insert(boundaryPositions.end(), pos.begin(), pos.end());
    }
    boundaryParticles.update(boundaryPositions, sf::Color::White);
    
    sim.addBoundary(boundaries);
    */

    enum State currentState = FLUID;
    enum BoundaryState currentBoundaryType = POLYGON;
    std::vector<Boundary> boundaries;
    std::vector<Vector2f> boundaryVertices;
    BoundaryData data = {
        boundaryVertices,
        boundaryRadius,
        1.f
    };

    while (window.isOpen()) {
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) window.close();

            if (const auto* keyPress = event->getIf<sf::Event::KeyPressed>()) {
                if (keyPress->scancode == sf::Keyboard::Scan::Num1) {
                    currentState = FLUID;
                } else if (keyPress->scancode == sf::Keyboard::Scan::Num2) {
                    currentState = BOUNDARY;
                }
            }

            if (currentState == FLUID) {
                handleFluidState(event, sim);
            } else if (currentState == BOUNDARY) {
                handleBoundaryState(event, sim, boundaries, data);
            }
        }
        
		std::vector<Vector2f> pos = sim.getPosition();
		std::vector<Vector2f> vel = sim.getVelocity();
        const auto d = sim.getDensity();
        particles.update(pos, vel, 0.f, 100.f, ColorMap::viridis);

        std::vector<Vector2f> boundaryPositions;
        for (const auto b : boundaries) {
            const auto pos = b.getBoundaryParticlePositions();
            boundaryPositions.insert(boundaryPositions.end(), pos.begin(), pos.end());
        }
        boundaryParticles.update(boundaryPositions, sf::Color::White);
        
        window.clear();
        window.draw(particles);
        window.draw(boundaryParticles);
        window.display();
        sim.update();
        SimFrameMark;
    }
}

void handleFluidState(std::optional<sf::Event> event, Simulation& sim) {
    if (const auto* mouseClick = event->getIf<sf::Event::MouseButtonPressed>()) {
        if (mouseClick->button == sf::Mouse::Button::Left) {
            printf("fleft\n");
        } else if (mouseClick->button == sf::Mouse::Button::Right) {
            printf("fright\n");
        }
    } else if (const auto* scroll = event->getIf<sf::Event::MouseWheelScrolled>()) {
        if (scroll->delta > 0.f) {
            printf("fscroll up\n");
        } else if (scroll->delta < 0.f) {
            printf("fscroll down\n");
        }
    }
}

void handleBoundaryState(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
) {
    if (const auto* mouseClick = event->getIf<sf::Event::MouseButtonPressed>()) {
        if (mouseClick->button == sf::Mouse::Button::Left) {
            // To stop the current vertex from being modified further, push a new one
            const Vector2f newVertPos = {
                (float) mouseClick->position.x, 
                (float) mouseClick->position.y
            };

            // If this is the first vertex, we need to save it first
            if (data.boundaryVertices.empty()) data.boundaryVertices.push_back(newVertPos);

            // Add a temporary vertex to be modified, at the same position
            data.boundaryVertices.push_back(newVertPos);
        } else if (mouseClick->button == sf::Mouse::Button::Right) {
            // Remove the temporary vertex and boundary
            data.boundaryVertices.pop_back();
            boundaries.pop_back();

            // Create the boundary to be saved
            Boundary b(data.smoothingRadius);
            b.createPolygon(data.boundaryVertices, data.compression);
            
            // Activate and add boundary for simulation
            b.activateBoundary();
            sim.addBoundary(b);
            boundaries.push_back(b);
            
            // Append a new boundary to be modified instead of current
            boundaries.push_back(b);

            data.boundaryVertices.clear();
        }
    } else if (const auto* scroll = event->getIf<sf::Event::MouseWheelScrolled>()) {
        data.compression += scroll->delta / 25.f;
        data.compression = std::clamp(data.compression, 0.f, 1.f);

        if (boundaries.empty()) return;

        // Update the compression by creating a new boundary
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createPolygon(data.boundaryVertices, data.compression);
        boundaries.back() = tempBoundary;
    } else if (const auto* mouseMove = event->getIf<sf::Event::MouseMoved>()) {
        if (data.boundaryVertices.empty()) return;

        // Modify the last vertex position when mouse moves
        auto& lastVert = data.boundaryVertices.back();
        lastVert = {
            (float) mouseMove->position.x, 
            (float) mouseMove->position.y
        };
        
        // Create a temporary boundary to draw current vertices
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createPolygon(data.boundaryVertices, data.compression);

        if (boundaries.empty()) {
            boundaries.push_back(tempBoundary);
        } else {
            boundaries.back() = tempBoundary;
        }
    }
}