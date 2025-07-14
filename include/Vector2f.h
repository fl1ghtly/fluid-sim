#pragma once
#include <math.h>

struct Vector2f {
	float x;
	float y;

	Vector2f& operator+=(const Vector2f& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		return *this;	
	}

	Vector2f& operator-=(const Vector2f& rhs) {
		this->x -= rhs.x;
		this->y -= rhs.y;
		return *this;
	}

	Vector2f& operator*=(const float& n) {
		this->x *= n;
		this->y *= n;
		return *this;
	}

	Vector2f& operator/=(const float& n) {
		this->x /= n;
		this->y /= n;
		return *this;
	}

	Vector2f normalize() const;
	float magnitude() const;
	float dot(const Vector2f&) const;
};

inline Vector2f operator+(const Vector2f& lhs, const Vector2f& rhs) {
	Vector2f res = lhs;
	res += rhs;
	return res;
}

inline Vector2f operator-(const Vector2f& lhs, const Vector2f& rhs) {
	Vector2f res = lhs;
	res -= rhs;
	return res;
}

inline Vector2f operator*(const Vector2f& lhs, const float& n) {
	Vector2f res = lhs;
	res *= n;
	return res;
}

inline Vector2f operator*(const float& n, const Vector2f& rhs) {
	Vector2f res = rhs;
	res *= n;
	return res;
}

inline Vector2f operator/(const Vector2f& lhs, const float& n) {
	Vector2f res = lhs;
	res /= n;
	return res;
}