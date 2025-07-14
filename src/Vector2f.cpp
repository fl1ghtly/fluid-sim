#include "Vector2f.h"

float Vector2f::magnitude() const {
	return sqrt(this->x * this->x + this->y * this->y);
}

Vector2f Vector2f::normalize() const{
	return *this / this->magnitude();
}

float Vector2f::dot(const Vector2f& rhs) const {
	return this->x * rhs.x + this->y * rhs.y;
}