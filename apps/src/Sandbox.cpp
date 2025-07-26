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

constexpr float SCROLL_TICK_RESOLUTION = 1.f / 25.f;

enum State {
    BOUNDARY,
    FLUID
};

enum FluidState {
    VELOCITY,
    DENSITY,
    PRESSURE
};

enum BoundaryState {
    POLYGON,
    BOX,
    CIRCLE
};

struct BoundaryData {
    Boundary& tempBoundary;
    std::vector<Vector2f>& boundaryVertices;
    float smoothingRadius;
    float compression;
    bool creatingBoundary = false;
};

struct FluidPlaceData {
    Vector2f origin;
    float radius;
};

void handleFluidState(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    FluidPlaceData& data
);
void handleBoundaryState(
    std::optional<sf::Event> event, 
    BoundaryState state,
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
);
void handlePolygon(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
);
void handleBox(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
);
void handleCircle(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
);
void stopDrawingTemporaryBoundary(BoundaryData& data);

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
    
	Simulation sim(width, height, params, downDir);
    ParticleSystem particles;
    ParticleSystem boundaryParticles;

    enum State currentState = FLUID;
    enum FluidState currentFluidState = VELOCITY;
    enum BoundaryState currentBoundaryType = POLYGON;

    FluidPlaceData fData = {
        {0.f, 0.f},
        100.f,
    };
    
    std::vector<Boundary> boundaries;
    std::vector<Vector2f> boundaryVertices;
    Boundary temp(boundaryRadius);
    BoundaryData bData = {
        temp,
        boundaryVertices,
        boundaryRadius,
        1.f,
    };

    while (window.isOpen()) {
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) window.close();

            if (const auto* keyPress = event->getIf<sf::Event::KeyPressed>()) {
                switch (keyPress->code) {
                    case sf::Keyboard::Key::Q:
                        currentState = FLUID;
                        stopDrawingTemporaryBoundary(bData);
                        break;
                    case sf::Keyboard::Key::E:
                        currentState = BOUNDARY;
                        stopDrawingTemporaryBoundary(bData);
                        break;
                    case sf::Keyboard::Key::Z:
                        currentBoundaryType = POLYGON;
                        stopDrawingTemporaryBoundary(bData);
                        break;
                    case sf::Keyboard::Key::X:
                        currentBoundaryType = BOX;
                        stopDrawingTemporaryBoundary(bData);
                        break;
                    case sf::Keyboard::Key::C:
                        currentBoundaryType = CIRCLE;
                        stopDrawingTemporaryBoundary(bData);
                        break;
                    case sf::Keyboard::Key::Num1:
                        currentFluidState = VELOCITY;
                        break;
                    case sf::Keyboard::Key::Num2:
                        currentFluidState = DENSITY;
                        break;
                    case sf::Keyboard::Key::Num3:
                        currentFluidState = PRESSURE;
                        break;
                    default:
                        break;
                }
            }

            if (currentState == FLUID) {
                handleFluidState(event, sim, fData);
            } else if (currentState == BOUNDARY) {
                handleBoundaryState(event, currentBoundaryType, sim, boundaries, bData);
            }
        }
        
        // Render fluid particles
		std::vector<Vector2f> pos = sim.getPosition();
        switch (currentFluidState) {
            case VELOCITY:
            {
                const auto vel = sim.getVelocity();
                particles.update(pos, vel, 0.f, 100.f, ColorMap::viridis);
                break;
            }
            case DENSITY:
            {
                const auto d = sim.getDensity();
                particles.update(pos, d, 0.f, 5E+3f, ColorMap::viridis);
                break;
            }
            case PRESSURE:
            {
                const auto p = sim.getPressure();
                particles.update(pos, p, 0.f, 1E+6f, ColorMap::viridis);
                break;
            }
        }

        // Render all boundaries
        std::vector<Vector2f> boundaryPositions;
        for (const auto b : boundaries) {
            const auto pos = b.getBoundaryParticlePositions();
            boundaryPositions.insert(boundaryPositions.end(), pos.begin(), pos.end());
        }
        const auto tempPos = bData.tempBoundary.getBoundaryParticlePositions();
        boundaryPositions.insert(boundaryPositions.end(), tempPos.begin(), tempPos.end());
        boundaryParticles.update(boundaryPositions, sf::Color::White);

        window.clear();
        window.draw(particles);
        window.draw(boundaryParticles);
        // Draw reticle to show where fluid will be placed
        if (currentState == FLUID) {
            sf::CircleShape c(fData.radius);
            c.setFillColor(sf::Color::Transparent);
            c.setOutlineColor(sf::Color::White);
            c.setOutlineThickness(0.5f);
            c.setOrigin({fData.radius, fData.radius});
            c.setPosition({fData.origin.x, fData.origin.y});
            window.draw(c);
        }
        window.display();
        sim.update();
        SimFrameMark;
    }
}

void handleFluidState(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    FluidPlaceData& data
) {
    if (const auto* mouseClick = event->getIf<sf::Event::MouseButtonPressed>()) {
        if (mouseClick->button == sf::Mouse::Button::Left) {
            sim.initializeParticleCircle(data.origin, data.radius);
        } else if (mouseClick->button == sf::Mouse::Button::Right) {
            printf("fright\n");
        }
    } else if (const auto* scroll = event->getIf<sf::Event::MouseWheelScrolled>()) {
        data.radius += scroll->delta;
        // Keep radius >= 0
        data.radius = std::max(0.f, data.radius);
    } else if (const auto* mouseMove = event->getIf<sf::Event::MouseMoved>()) {
        // Change where fluid will be affected
        data.origin = {
            (float) mouseMove->position.x,
            (float) mouseMove->position.y
        };
    }
}

void handleBoundaryState(
    std::optional<sf::Event> event, 
    BoundaryState state,
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
) {
    switch (state)
    {
        case POLYGON:
            handlePolygon(event, sim, boundaries, data);
            break;
        case BOX:
            handleBox(event, sim, boundaries, data);
            break;
        case CIRCLE:
            handleCircle(event, sim, boundaries, data);
            break;
    }
}

void handlePolygon(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
) {
    if (const auto* mouseClick = event->getIf<sf::Event::MouseButtonPressed>()) {
        if (mouseClick->button == sf::Mouse::Button::Left) {
            // Save current mouse position
            const Vector2f newVertPos = {
                (float) mouseClick->position.x, 
                (float) mouseClick->position.y
            };
            
            // If this is the first vertex, we need to save it first
            if (!data.creatingBoundary) {
                data.boundaryVertices.push_back(newVertPos);
                data.creatingBoundary = true;
            }
    
            // Add a temporary vertex to be modified, at the same position
            data.boundaryVertices.push_back(newVertPos);
        } else if (mouseClick->button == sf::Mouse::Button::Right && data.creatingBoundary) {
            // Remove the temporary vertex and boundary
            data.boundaryVertices.pop_back();
            data.tempBoundary = Boundary(data.smoothingRadius);
    
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
            data.creatingBoundary = false;
        }
    } else if (const auto* scroll = event->getIf<sf::Event::MouseWheelScrolled>()) {
        data.compression += scroll->delta * SCROLL_TICK_RESOLUTION;
        data.compression = std::clamp(data.compression, 0.f, 1.f);
    
        if (!data.creatingBoundary) return;
    
        // Update the compression by creating a new boundary
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createPolygon(data.boundaryVertices, data.compression);
        data.tempBoundary = tempBoundary;
    } else if (const auto* mouseMove = event->getIf<sf::Event::MouseMoved>()) {
        if (!data.creatingBoundary) return;
    
        // Modify the last vertex position when mouse moves
        auto& lastVert = data.boundaryVertices.back();
        lastVert = {
            (float) mouseMove->position.x, 
            (float) mouseMove->position.y
        };
        
        // Create a temporary boundary to draw current vertices
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createPolygon(data.boundaryVertices, data.compression);
        data.tempBoundary = tempBoundary;
    }
}

void handleBox(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
) {
    if (const auto* mouseClick = event->getIf<sf::Event::MouseButtonPressed>()) {
        if (mouseClick->button == sf::Mouse::Button::Left) {
            // Save vertex position on left click
            const Vector2f newVertPos = {
                (float) mouseClick->position.x, 
                (float) mouseClick->position.y
            };
            data.boundaryVertices.push_back(newVertPos);
            
            // Check if creating a new boundary
            if (!data.creatingBoundary) {
                data.creatingBoundary = true;
                return;
            }
            
            // Already creating boundary, therefore this is the second vertex
            // and we can create the box
            Boundary b(data.smoothingRadius);
            b.createBox(
                data.boundaryVertices[0], 
                data.boundaryVertices[1], 
                data.compression
            );

            // Remove temp boundary
            data.tempBoundary = Boundary(data.smoothingRadius);

            b.activateBoundary();
            sim.addBoundary(b);
            boundaries.push_back(b);
            data.boundaryVertices.clear();
            data.creatingBoundary = false;
        } else if (mouseClick->button == sf::Mouse::Button::Right) {
            // Cancel creation of box boundary
            data.boundaryVertices.clear();
            data.tempBoundary = Boundary(data.smoothingRadius);
            data.creatingBoundary = false;
        }
    } else if (const auto* scroll = event->getIf<sf::Event::MouseWheelScrolled>()) {
        data.compression += scroll->delta * SCROLL_TICK_RESOLUTION;
        data.compression = std::clamp(data.compression, 0.f, 1.f);
    
        if (!data.creatingBoundary) return;

        const Vector2f lastVert = {
            (float) scroll->position.x, 
            (float) scroll->position.y
        };

        // Update the compression by creating a new boundary
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createBox(
            data.boundaryVertices[0], 
            lastVert,
            data.compression
        );
        data.tempBoundary = tempBoundary;
    } else if (const auto* mouseMove = event->getIf<sf::Event::MouseMoved>()) {
        if (!data.creatingBoundary) return;
    
        // Modify the last vertex position when mouse moves
        const Vector2f lastVert = {
            (float) mouseMove->position.x, 
            (float) mouseMove->position.y
        };
        
        // Create a temporary boundary to draw current vertices
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createBox(
            data.boundaryVertices[0], 
            lastVert, 
            data.compression
        );
        data.tempBoundary = tempBoundary;
    }
}

void handleCircle(
    std::optional<sf::Event> event, 
    Simulation& sim, 
    std::vector<Boundary>& boundaries,
    BoundaryData& data
) {
    if (const auto* mouseClick = event->getIf<sf::Event::MouseButtonPressed>()) {
        if (mouseClick->button == sf::Mouse::Button::Left) {
            const Vector2f newVertPos = {
                (float) mouseClick->position.x, 
                (float) mouseClick->position.y
            };
            
            // Check if creating a new boundary
            if (!data.creatingBoundary) {
                // Save only 1 vertex position for circle
                data.boundaryVertices.push_back(newVertPos);
                data.creatingBoundary = true;
                return;
            }
    
            // Already creating boundary, therefore we finish and save circle boundary
            Boundary b(data.smoothingRadius);
            b.createCircle(
                data.boundaryVertices[0], 
                (newVertPos - data.boundaryVertices[0]).magnitude(), 
                data.compression
            );
            
            // Remove temporary boundary
            data.tempBoundary = Boundary(data.smoothingRadius);

            b.activateBoundary();
            sim.addBoundary(b);
            boundaries.push_back(b);
            data.boundaryVertices.clear();
            data.creatingBoundary = false;
        } else if (mouseClick->button == sf::Mouse::Button::Right) {
            // Cancel creation of boundary
            data.boundaryVertices.clear();
            data.tempBoundary = Boundary(data.smoothingRadius);
            data.creatingBoundary = false;
        }
    } else if (const auto* scroll = event->getIf<sf::Event::MouseWheelScrolled>()) {
        data.compression += scroll->delta * SCROLL_TICK_RESOLUTION;
        data.compression = std::clamp(data.compression, 0.f, 1.f);
    
        if (!data.creatingBoundary) return;
    
        const Vector2f lastVert = {
            (float) scroll->position.x, 
            (float) scroll->position.y
        };
    
        // Update the compression by creating a new boundary
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createCircle(
            data.boundaryVertices[0], 
            (lastVert - data.boundaryVertices[0]).magnitude(),
            data.compression
        );
        data.tempBoundary = tempBoundary;
    } else if (const auto* mouseMove = event->getIf<sf::Event::MouseMoved>()) {
        if (!data.creatingBoundary) return;
    
        // Modify the last vertex position when mouse moves
        const Vector2f lastVert = {
            (float) mouseMove->position.x, 
            (float) mouseMove->position.y
        };
        
        // Create a temporary boundary to draw current vertices
        Boundary tempBoundary(data.smoothingRadius);
        tempBoundary.createCircle(
            data.boundaryVertices[0], 
            (lastVert - data.boundaryVertices[0]).magnitude(), 
            data.compression
        );
        data.tempBoundary = tempBoundary;
    }
}

void stopDrawingTemporaryBoundary(BoundaryData& data) {
    if (data.creatingBoundary) {
        data.boundaryVertices.clear();
        data.tempBoundary = Boundary(data.smoothingRadius);
        data.creatingBoundary = false;
    } 
}