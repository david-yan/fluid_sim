#include <math.h>

#include "vector2f.hpp"

float Vector2f::magnitude()
{
    return std::sqrt(x * x + y * y);
}

Vector2f Vector2f::normalized()
{
    float mag = magnitude();
    return Vector2f(x / mag, y / mag);
}

float Vector2f::inner_product(const Vector2f &other)
{
    return x * other.x + y * other.y;
}

Vector2f randomUnitVector2f() {
    float theta = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * 2.f * M_PI;
    return Vector2f(sin(theta), cos(theta));
}