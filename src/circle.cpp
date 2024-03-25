#include <math.h>
#include <SFML/Graphics/CircleShape.hpp>

#include "circle.hpp"
#include "vector2f.hpp"

Vector2f Circle::getCenter()
{
    return getPosition() + centerDisplacement;
}

void Circle::setCenter(Vector2f center)
{
    setPosition(center - centerDisplacement);
}

float Circle::distanceFromCenter(Vector2f point)
{
    Vector2f diff = getCenter() - point;
    return diff.magnitude();
}