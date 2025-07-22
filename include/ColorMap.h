#pragma once
#include <vector>
#include <SFML/Graphics.hpp>

namespace ColorMap {
	const std::vector<sf::Color> viridis = {
		{68, 1, 84},
		{65, 68, 135},
		{42, 120, 142},
		{34, 168, 132},
		{122, 209, 81},
		{253, 231, 37},
	};

	const std::vector<sf::Color> inferno = {
		{0, 0, 4},
		{66, 10, 104},
		{147, 38, 103},
		{221, 81, 58},
		{252, 165, 10},
		{252, 255, 164},
	};

	const std::vector<sf::Color> magma = {
		{0, 0, 4},
		{59, 15, 112},
		{140, 41, 129},
		{222, 73, 104},
		{254, 159, 109},
		{252, 253, 191},
	};
}