#pragma once

#include <SFML/Graphics.hpp>

class Vector2f: public sf::Vector2f {
public:
    Vector2f(): sf::Vector2f() {};
    Vector2f(float x, float y): sf::Vector2f(x, y) {};
    Vector2f(sf::Vector2f &vec): sf::Vector2f(vec.x, vec.y) {};

    float magnitude();
    Vector2f normalized();
    float inner_product(Vector2f &other);
};