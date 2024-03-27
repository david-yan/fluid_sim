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
    std::vector<unsigned int> _getNeighboringQuadrants(const unsigned int quadrant);

public:
    std::vector<Circle> particles;
    std::vector<Vector2f> particleVelocities;
    std::vector<float> particleDensities;
    std::vector<Vector2f> particleForces;
    
    std::vector<std::vector<float>> allDensities;
    std::vector<std::vector<float>> allPressures;

    ParticleStore();
    ParticleStore(sf::Vector2u &bounds, const float smoothingRadius)
    : smoothingRadius(smoothingRadius) {
        width = bounds.x;
        height = bounds.y;
        allDensities.resize(width / DENSITY_SCALE, std::vector<float>(height / DENSITY_SCALE));
        allPressures.resize(width / DENSITY_SCALE, std::vector<float>(height / DENSITY_SCALE));
    };

    void spawnParticles(const int rows, const int cols, const float radius, const Vector2f &center, const float spacing = 0.f);
    void spawnRandomParticles(const int numParticles, const float radius, const sf::Vector2u bounds);

    void applyForces(sf::Time &elapsed);

    void setSmoothingRadius(const float newSmoothingRadius);
    void updateQuadrants();
    void visualizeQuadrants(sf::RenderWindow &window);

    float calculateDensity(const Vector2f &point);
    float calculatePressure(const float density);

    void calculateAllDensities();
    void calculateAllParticleDensities();
    void visualizeDensity(sf::RenderWindow &window);

    void calculateAllPressures();
    void visualizePressure(sf::RenderWindow &window);

    Vector2f calculateParticleForce(int idx);
    void calculateParticleForces();
    void visualizeParticleForces(sf::RenderWindow &window);

    void drawParticles(sf::RenderWindow &window);
};