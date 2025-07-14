#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "fluid.h"
#include "Vector2f.h"

int main(void) {
	auto window = sf::RenderWindow(sf::VideoMode({1920u, 1080u}), "CMake SFML Project");
    window.setFramerateLimit(144);
	
	Fluid fluid(1920, 1080, 100, 5);
    while (window.isOpen())
    {
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }

        window.clear();
        window.display();
    }
}