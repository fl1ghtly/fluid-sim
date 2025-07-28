#include "GridCell.h"

GridCell::GridCell(float x, float y, float cellSize, Vector2f min): cellSize(cellSize) {
	this->x = static_cast<std::size_t>((x - min.x) / cellSize);
	this->y = static_cast<std::size_t>((y - min.y) / cellSize);
	calculateZOrder();
}

GridCell::GridCell(int x, int y, float cellSize, Vector2f min): 
	cellSize(cellSize) 
{
	const float xF = static_cast<float>(x);
	const float yF = static_cast<float>(y);
	this->x = static_cast<std::size_t>(xF - min.x);
	this->y = static_cast<std::size_t>(yF - min.y);
	calculateZOrder();
}

GridCell::GridCell(Vector2f pos, float cellSize, Vector2f min): cellSize(cellSize) {
	this->x = static_cast<std::size_t>((pos.x - min.x) / cellSize);
	this->y = static_cast<std::size_t>((pos.y - min.y) / cellSize);
	calculateZOrder();
}

bool GridCell::operator==(const GridCell& other) const {
	return x == other.x && y == other.y;
}

void GridCell::calculateZOrder() {
	constexpr std::size_t bitmask[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
	constexpr std::size_t shift[] = {1, 2, 4, 8};

	// Interleave bits from Y and X, starting with Y to get the Morton Code
	std::size_t uX = (x | (x << shift[3u])) & bitmask[3u];
	uX = (uX | (uX << shift[2])) & bitmask[2];
	uX = (uX | (uX << shift[1])) & bitmask[1];
	uX = (uX | (uX << shift[0])) & bitmask[0];

	std::size_t uY = (y | (y << shift[3u])) & bitmask[3u];
	uY = (uY | (uY << shift[2])) & bitmask[2];
	uY = (uY | (uY << shift[1])) & bitmask[1];
	uY = (uY | (uY << shift[0])) & bitmask[0];

	zOrder = uX | (uY << 1);
}

std::size_t std::hash<GridCell>::operator()(const GridCell& cell) const {
	// Prime Numbers
	constexpr std::size_t p1 = 73856093;
	constexpr std::size_t p2 = 19349663;

	return (cell.x * p1) ^ (cell.y * p2);
}
