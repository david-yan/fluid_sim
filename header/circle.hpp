#pragma once

#include <SFML/Graphics/CircleShape.hpp>

#include "vector2f.hpp"

class Circle: public sf::CircleShape {
protected:
    Vector2f centerDisplacement;
    
public:
    Circle(float radius, std::size_t pointCount=30): sf::CircleShape(radius, pointCount) {
        centerDisplacement.x = radius;
        centerDisplacement.y = radius;
    };
    Circle(Vector2f center, float radius, std::size_t pointCount=30): sf::CircleShape(radius, pointCount) {
        centerDisplacement.x = radius;
        centerDisplacement.y = radius;
        setPosition(center - centerDisplacement);
    };

    Vector2f getCenter();
    void setCenter(Vector2f &center);

    float distanceFromCenter(const Vector2f &point);
};