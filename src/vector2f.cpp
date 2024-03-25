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

float Vector2f::inner_product(Vector2f &other)
{
    return x * other.x + y * other.y;
}
