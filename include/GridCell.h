#pragma once
#include <unordered_map>
#include <Vector2f.h>

struct GridCell {
	public:
		int x, y;
		unsigned int zOrder;
		float cellSize;

		GridCell(float x, float y, float cellSize);
		GridCell(int x, int y, float cellSize);
		GridCell(Vector2f pos, float cellSize);
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