#include <iostream>

#include "particles.hpp"
#include "circle.hpp"
#include "kernels.hpp"
#include "timer.hpp"
#include "vector2f.hpp"

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
    updateQuadrants();
}

void ParticleStore::spawnRandomParticles(const int numParticles, const float radius, const sf::Vector2u bounds) {
    std::srand(std::time(nullptr));
    for (int i = 0; i < numParticles; i++) {
        float x = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * static_cast<float>(bounds.x - (2 * radius));
        float y = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * static_cast<float>(bounds.y - (2 * radius));

        Circle particle(Vector2f(x, y), radius);
        particle.setFillColor(sf::Color(255, 255, 255));
        particles.push_back(particle);
        particleVelocities.push_back(Vector2f(0.f, 0.f));
    }
    updateQuadrants();
}

void ParticleStore::applyForces(sf::Time &elapsed) {
    for (int i = 0; i < particleForces.size(); i++) {
        Circle &particle = particles[i];
        Vector2f &particleVelocity = particleVelocities[i];

        float radius = particle.getRadius();

        // Update velocity
        // particleVelocity += GRAVITY * elapsed.asSeconds();
        Vector2f acc = particleForces[i] / particleDensities[i];
        particleVelocity += acc * elapsed.asSeconds();

        if (particleVelocity.magnitude() > MAX_VELOCITY) {
            particleVelocity *= MAX_VELOCITY / particleVelocity.magnitude();
        }

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
    }
    updateQuadrants();
}

void ParticleStore::setSmoothingRadius(const float newSmoothingRadius) {
    smoothingRadius = newSmoothingRadius;
    updateQuadrants();
}

unsigned int ParticleStore::_positionToQuadrant(const Vector2f &pos) {
    return static_cast<unsigned int>(pos.y / smoothingRadius) * (smoothingRadius + 1) + static_cast<unsigned int>(pos.x / smoothingRadius);
}

unsigned int ParticleStore::_getParticleQuadrant(const unsigned int idx){
    const Vector2f &center = particles[idx].getCenter();
    return _positionToQuadrant(center);
}

std::vector<unsigned int> ParticleStore::_getNeighboringQuadrants(const unsigned int quadrant) {
    unsigned int xLimit = (width / smoothingRadius) + 1;
    unsigned int yLimit = (height / smoothingRadius) + 1;

    unsigned int x = quadrant % (static_cast<int>(smoothingRadius) + 1);
    unsigned int y = quadrant / (static_cast<int>(smoothingRadius) + 1);

    std::vector<unsigned int> neighbors;
    int deltas[] = {-1, 0, 1};
    for (int dx : deltas) {
        for (int dy : deltas) {
            unsigned int newX = x + dx;
            unsigned int newY = y + dy;

            if (newX < 0 || newY < 0 || newX >= xLimit || newY >= yLimit) {
                continue;
            }
            neighbors.push_back(newY * (smoothingRadius + 1) + newX);
        }
    }
    return neighbors;
}

void ParticleStore::updateQuadrants() {
    quadrantToParticleIdxs.clear();

    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        unsigned int quadrant = _getParticleQuadrant(i);
        quadrantToParticleIdxs[quadrant].push_back(i);
    }
}

void ParticleStore::visualizeQuadrants(sf::RenderWindow &window) {
    for (int i = 0; i < width / smoothingRadius; i++) {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f((i+1)*smoothingRadius, 0.f), sf::Color(57, 57, 57)),
            sf::Vertex(sf::Vector2f((i+1)*smoothingRadius, height), sf::Color(57, 57, 57))
        };
        window.draw(line, 2, sf::LineStrip);
    }
    for (int j = 0; j < height / smoothingRadius; j++) {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(0.f, (j+1)*smoothingRadius), sf::Color(57, 57, 57)),
            sf::Vertex(sf::Vector2f(width, (j+1)*smoothingRadius), sf::Color(57, 57, 57))
        };
        window.draw(line, 2, sf::LineStrip);
    }
}

float ParticleStore::calculateDensity(const Vector2f &point) {
    float density = 0;

    unsigned int quadrant = _positionToQuadrant(point);
    std::vector<unsigned int> neighboringQuadrants = _getNeighboringQuadrants(quadrant);

    for (int quadrant : neighboringQuadrants) {
        std::vector<unsigned int> &particleIdxs = quadrantToParticleIdxs[quadrant];
        for (int i : particleIdxs) {
            Circle &particle = particles[i];
            const float dst = particle.distanceFromCenter(point);
            // Visualize the particles in neighboring quadrants
            // if (dst < smoothingRadius) {
            //     particle.setFillColor(sf::Color(0, 0, 255));
            // }
            // else {
            //     particle.setFillColor(sf::Color(255, 0, 0));
            // }

            const float influence = spikyPow2Kernel(smoothingRadius, dst);
            density += MASS * influence;
        }
    }
    return density;
}

float ParticleStore::calculatePressure(const float density) {
    // std::cout << "density: " << density << std::endl;
    float densityError = density - TARGET_DENSITY;
    // std::cout << "densityError: " << densityError << std::endl;
    float pressure = PRESSURE_MULTIPLIER * densityError;
    return pressure;
}

void ParticleStore::calculateAllDensities() {
    int i, j;
    float density;
    
    // Timer timer;

    #pragma omp parallel for default(none) \
                             private(i, j, density) shared(allDensities, DENSITY_SCALE)
    for (i = 0; i < allDensities.size(); i++) {
        for (j = 0; j < allDensities[0].size(); j++) {
            const Vector2f point(i * DENSITY_SCALE, j * DENSITY_SCALE);
            float density = calculateDensity(point);
            allDensities[i][j] = density;
        }
    }

    // std::cout << "Time: " << timer.elapsed() << std::endl;
}

void ParticleStore::calculateAllParticleDensities() {
    particleDensities.clear();
    
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        float density = calculateDensity(particles[i].getCenter());
        particleDensities.push_back(density);
    }
}

Vector2f ParticleStore::calculateParticleForce(int idx) {
    Vector2f gradient;

    Circle &particle = particles[idx];
    const Vector2f center = particle.getCenter();
    unsigned int quadrant = _positionToQuadrant(center);
    std::vector<unsigned int> neighboringQuadrants = _getNeighboringQuadrants(quadrant);

    for (int quadrant : neighboringQuadrants) {
        std::vector<unsigned int> &particleIdxs = quadrantToParticleIdxs[quadrant];
        for (int i : particleIdxs) {
            if (i == idx) continue;

            float dst = Vector2f(center - particles[i].getCenter()).magnitude();

            Vector2f dir = (dst != 0) ? (particles[i].getCenter() - center) / dst : randomUnitVector2f();
            float slope = spikyPow2KernelDerivative(smoothingRadius, dst);
            float density = particleDensities[i];

            float sharedPressure = (calculatePressure(density) + calculatePressure(particleDensities[idx])) / 2.f;
            gradient += calculatePressure(density) * dir * slope * 100.f / density;
        }
    }
    return gradient;
}

void ParticleStore::calculateParticleForces() {
    particleForces.clear();

    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        Vector2f pressure = calculateParticleForce(i);

        // TODO: Fix this scaling issue
        particleForces.push_back(pressure);
        // particleForces.push_back(pressure.normalized() * 20.f);
    }
}

void ParticleStore::visualizeParticleForces(sf::RenderWindow &window) {
    #pragma omp parallel for
    for (int i = 0; i < particleForces.size(); i++) {
        const Vector2f &pos = particles[i].getCenter();
        const Vector2f &force = particleForces[i];

        sf::Vertex line[] = {
            sf::Vertex(pos, sf::Color(255, 0, 0)),
            sf::Vertex(pos + force, sf::Color(255, 0, 0))
        };
        window.draw(line, 2, sf::LineStrip);
    }
}

void ParticleStore::visualizeDensity(sf::RenderWindow &window) {
    #pragma omp parallel for
    for (int i = 0; i < allDensities.size(); i++) {
        for (int j = 0; j < allDensities[0].size(); j++) {
            const float density = allDensities[i][j];
            sf::RectangleShape pixel(Vector2f(DENSITY_SCALE, DENSITY_SCALE));
            pixel.setPosition(Vector2f(i*DENSITY_SCALE, j*DENSITY_SCALE));
            pixel.setFillColor(sf::Color(0, 0, 255, 255 * std::min(density / 2, 1.f)));
            window.draw(pixel);
        }
    }
}

void ParticleStore::calculateAllPressures() {
    int i, j;
    float density;
    
    // Timer timer;

    #pragma omp parallel for default(none) \
                             private(i, j, density, pressure) shared(allPressures, DENSITY_SCALE)
    for (i = 0; i < allPressures.size(); i++) {
        for (j = 0; j < allPressures[0].size(); j++) {
            const Vector2f point(i * DENSITY_SCALE, j * DENSITY_SCALE);
            float density = calculateDensity(point);
            float pressure = calculatePressure(density);
            // std::cout << "pressure: " << pressure << std::endl;
            allPressures[i][j] = pressure;
        }
    }
}

void ParticleStore::visualizePressure(sf::RenderWindow &window){
    #pragma omp parallel for
    for (int i = 0; i < allPressures.size(); i++) {
        for (int j = 0; j < allPressures[0].size(); j++) {
            const float pressure = allPressures[i][j];
            // std::cout << "pressure: " << pressure << std::endl;
            float r = (pressure > 0) ? 0 : 255 * std::min(-pressure * 2, 1.f);
            float b = (pressure < 0) ? 0 : 255 * std::min(pressure * 2, 1.f);

            sf::RectangleShape pixel(Vector2f(DENSITY_SCALE, DENSITY_SCALE));
            pixel.setPosition(Vector2f(i*DENSITY_SCALE, j*DENSITY_SCALE));
            pixel.setFillColor(sf::Color(r, 0, b, 255));
            window.draw(pixel);
        }
    }
}

void ParticleStore::drawParticles(sf::RenderWindow &window) {
    for (int i = 0; i < particles.size(); i++) {
        Circle &particle = particles[i];
        window.draw(particle);
    }
}