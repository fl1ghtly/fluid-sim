#include "GridCell.h"

GridCell::GridCell(float x, float y, float cellSize): cellSize(cellSize) {
	this->x = std::floor(x / cellSize);
	this->y = std::floor(y / cellSize);
}

GridCell::GridCell(int x, int y, float cellSize): 
	x(x),
	y(y),
	cellSize(cellSize) {}

GridCell::GridCell(Vector2f pos, float cellSize): cellSize(cellSize) {
	this->x = std::floor(pos.x / cellSize);
	this->y = std::floor(pos.y / cellSize);
}

bool GridCell::operator==(const GridCell& other) const {
	return x == other.x && y == other.y;
}

std::size_t std::hash<GridCell>::operator()(const GridCell& cell) const {
	// Prime Numbers
	constexpr int p1 = 73856093;
	constexpr int p2 = 19349663;

	return (cell.x * p1) ^ (cell.y * p2);
}
