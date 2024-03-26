#include <iostream>

#include "particles.hpp"
#include "circle.hpp"
#include "kernels.hpp"
#include "timer.hpp"

void ParticleStore::spawnParticles(const int rows, const int cols, const float radius, const Vector2f &center, const float spacing) {
    int nParticles = rows * cols;
    int rowCenterIdx = rows / 2;
    int colCenterIdx = cols / 2;
    if (rows % 2 == 1) {
        rowCenterIdx += 1;
    }
    if (cols % 2 == 1) {
        colCenterIdx += 1;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            Vector2f displacement((i - rowCenterIdx) * (2 * radius + spacing), (j - colCenterIdx) * (2 * radius + spacing));
            Circle particle(center + displacement, radius);
            particle.setFillColor(sf::Color(255, 255, 255));

            particles.push_back(particle);
            particleVelocities.push_back(Vector2f(0.f, 0.f));
        }
    }
}

void ParticleStore::spawnRandomParticles(const int numParticles, const float radius, const sf::Vector2u bounds) {
    std::srand(std::time(nullptr));
    for (int i = 0; i < numParticles; i++) {
        float x = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX) * static_cast <float> (bounds.x - (2 * radius));
        float y = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX) * static_cast <float> (bounds.y - (2 * radius));

        Circle particle(Vector2f(x, y), radius);
        particle.setFillColor(sf::Color(255, 255, 255));
        particles.push_back(particle);
        particleVelocities.push_back(Vector2f(0.f, 0.f));
    }
}

void ParticleStore::applyGravity(Circle &particle, Vector2f &particleVelocity, sf::Time &elapsed) {
    float radius = particle.getRadius();

    particle.move(particleVelocity);

    // Check particle bounds
    Vector2f updatedPosition = particle.getPosition();
    if (updatedPosition.x < 0 || updatedPosition.x + (2 * radius) >= width) {
        particleVelocity.x *= -1 * DAMPING.x;
    }
    if (updatedPosition.y < 0 || updatedPosition.y + (2 * radius) >= height) {
        particleVelocity.y *= -1 * DAMPING.y;
    }
    float boundedX = std::min(std::max(updatedPosition.x, 0.f), static_cast<float>(width - (2 * radius)));
    float boundedY = std::min(std::max(updatedPosition.y, 0.f), static_cast<float>(height - (2 * radius)));
    particle.setPosition(boundedX, boundedY);

    particleVelocity += GRAVITY * elapsed.asSeconds();
}

void ParticleStore::setSmoothingRadius(const float newSmoothingRadius) {
    smoothingRadius = newSmoothingRadius;
}

unsigned int ParticleStore::_positionToQuadrant(const Vector2f &pos) {
    return static_cast<unsigned int>(pos.y / smoothingRadius) * (smoothingRadius + 1) + static_cast<unsigned int>(pos.x / smoothingRadius);
}

unsigned int ParticleStore::_getParticleQuadrant(const unsigned int idx){
    const Vector2f &center = particles[idx].getCenter();
    return _positionToQuadrant(center);
}

void ParticleStore::updateQuadrants() {
    quadrantToParticleIdxs.clear();
    for (int i = 0; i < particles.size(); i++) {
        unsigned int quadrant = _getParticleQuadrant(i);
        if (quadrantToParticleIdxs.find(quadrant) != quadrantToParticleIdxs.end()) {
            quadrantToParticleIdxs[quadrant] = std::vector<unsigned int>(i);
        }
        else {
            quadrantToParticleIdxs[quadrant].push_back(i);
        }
    }
}

void ParticleStore::visualizeQuadrants(sf::RenderWindow &window) {
    
}

float ParticleStore::calculateDensity(const Vector2f &point, const float smoothingRadius) {
    const float mass = 1.f;
    float density = 0;

    for (int i = 0; i < particles.size(); i++) {
        Circle &particle = particles[i];
        const float dst = particle.distanceFromCenter(point);
        if (dst < smoothingRadius) {
            particle.setFillColor(sf::Color(0, 0, 255));
        }
        else {
            particle.setFillColor(sf::Color(255, 255, 255));
        }

        const float influence = smoothingKernel(smoothingRadius, dst);
        density += mass * influence;
    }
    return density;
}

void ParticleStore::calculateAllDensities(const float smoothingRadius) {
    int i, j;
    float density;
    
    Timer timer;

    #pragma omp parallel for default(none) \
                             private(i, j, density) shared(allDensities, DENSITY_SCALE)
    for (i = 0; i < allDensities.size(); i++) {
        for (j = 0; j < allDensities[0].size(); j++) {
            const Vector2f point(i * DENSITY_SCALE, j * DENSITY_SCALE);
            float density = calculateDensity(point, smoothingRadius);
            allDensities[i][j] = density;
        }
    }

    std::cout << "Time: " << timer.elapsed() << std::endl;
}

void ParticleStore::updateGravity(sf::Time &elapsed) {
    for (int i = 0; i < particles.size(); i++) {
        Circle &particle = particles[i];
        Vector2f &particleVelocity = particleVelocities[i];

        applyGravity(particle, particleVelocity, elapsed);
    }
}

void ParticleStore::drawParticles(sf::RenderWindow &window) {
    for (int i = 0; i < particles.size(); i++) {
        Circle &particle = particles[i];
        window.draw(particle);
    }
}