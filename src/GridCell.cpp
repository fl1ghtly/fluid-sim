#include "GridCell.h"

GridCell::GridCell(float x, float y, float cellSize): cellSize(cellSize) {
	this->x = std::floor(x / cellSize);
	this->y = std::floor(y / cellSize);
	calculateZOrder();
}

GridCell::GridCell(int x, int y, float cellSize): 
	x(x),
	y(y),
	cellSize(cellSize) 
{
	calculateZOrder();
}

GridCell::GridCell(Vector2f pos, float cellSize): cellSize(cellSize) {
	this->x = std::floor(pos.x / cellSize);
	this->y = std::floor(pos.y / cellSize);
	calculateZOrder();
}

bool GridCell::operator==(const GridCell& other) const {
	return x == other.x && y == other.y;
}

void GridCell::calculateZOrder() {
	constexpr unsigned int bitmask[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
	constexpr unsigned int shift[] = {1, 2, 4, 8};

	// Interleave bits from Y and X, starting with Y to get the Morton Code
	unsigned int uX = (x | (x << shift[3])) & bitmask[3];
	uX = (uX | (uX << shift[2])) & bitmask[2];
	uX = (uX | (uX << shift[1])) & bitmask[1];
	uX = (uX | (uX << shift[0])) & bitmask[0];

	unsigned int uY = (y | (y << shift[3])) & bitmask[3];
	uY = (uY | (uY << shift[2])) & bitmask[2];
	uY = (uY | (uY << shift[1])) & bitmask[1];
	uY = (uY | (uY << shift[0])) & bitmask[0];

	zOrder = uX | (uY << 1);
}

std::size_t std::hash<GridCell>::operator()(const GridCell& cell) const {
	// Prime Numbers
	constexpr int p1 = 73856093;
	constexpr int p2 = 19349663;

	return (cell.x * p1) ^ (cell.y * p2);
}
