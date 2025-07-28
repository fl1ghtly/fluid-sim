#pragma once
#include <unordered_map>
#include <Vector2f.h>

struct GridCell {
	public:
		std::size_t x, y;
		std::size_t zOrder;
		float cellSize;

		GridCell(float x, float y, float cellSize, Vector2f min={0.f, 0.f});
		GridCell(int x, int y, float cellSize, Vector2f min={0.f, 0.f});
		GridCell(Vector2f pos, float cellSize, Vector2f min={0.f, 0.f});
		bool operator==(const GridCell& other) const; 
	private:
		void calculateZOrder();
};

inline bool operator<(const GridCell& lhs, const GridCell& rhs) {
	return lhs.zOrder < rhs.zOrder;
}

template<>
struct std::hash<GridCell> {
	std::size_t operator()(const GridCell& cell) const;
};