#pragma once

#include <SFML/Graphics.hpp>
#include <unordered_map>

#include "constants.hpp"
#include "circle.hpp"
#include "vector2f.hpp"

class ParticleStore {
private:
    unsigned int width;
    unsigned int height;

    float smoothingRadius;
    std::unordered_map<unsigned int, std::vector<unsigned int>> quadrantToParticleIdxs;

    unsigned int _positionToQuadrant(const Vector2f &pos);
    unsigned int _getParticleQuadrant(const unsigned int idx);
public:
    std::vector<Circle> particles;
    std::vector<Vector2f> particleVelocities;

    std::vector<std::vector<float>> allDensities;

    ParticleStore();
    ParticleStore(sf::Vector2u &bounds, const float smoothingRadius)
    : smoothingRadius(smoothingRadius) {
        width = bounds.x;
        height = bounds.y;
        allDensities.resize(width / DENSITY_SCALE, std::vector<float>(height / DENSITY_SCALE));
    };

    void spawnParticles(const int rows, const int cols, const float radius, const Vector2f &center, const float spacing = 0.f);
    void spawnRandomParticles(const int numParticles, const float radius, const sf::Vector2u bounds);

    void applyGravity(Circle &particle, Vector2f &particleVelocity, sf::Time &elapsed);

    void setSmoothingRadius(const float newSmoothingRadius);
    void updateQuadrants();
    void visualizeQuadrants(sf::RenderWindow &window);

    float calculateDensity(const Vector2f &point, const float smoothingRadius);
    void calculateAllDensities(const float smoothingRadius);

    void drawParticles(sf::RenderWindow &window);
    void updateGravity(sf::Time &elapsed);
};