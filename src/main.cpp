#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>

#include "circle.hpp"
#include "vector2f.hpp"
#include "kernels.hpp"

Vector2f gravity(0.f, 9.8);
Vector2f damping(0.9, 0.9);

float minSmoothingRadius(1.0);
float smoothingRadius(50.0);

std::vector<Circle> particles;
std::vector<Vector2f> particleVelocities;

void spawnParticles(const int rows, const int cols, const float radius, const Vector2f &center, const float spacing = 0.f) {
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

void spawnRandomParticles(const int numParticles, const float radius, const sf::Vector2u bounds) {
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

void applyGravity(Circle &particle, Vector2f &particleVelocity, sf::Time &elapsed, const sf::Vector2u &bounds) {
    unsigned int width = bounds.x;
    unsigned int height = bounds.y;
    float radius = particle.getRadius();

    particle.move(particleVelocity);

    // Check particle bounds
    Vector2f updatedPosition = particle.getPosition();
    if (updatedPosition.x < 0 || updatedPosition.x + (2 * radius) >= width) {
        particleVelocity.x *= -1 * damping.x;
    }
    if (updatedPosition.y < 0 || updatedPosition.y + (2 * radius) >= height) {
        particleVelocity.y *= -1 * damping.y;
    }
    float boundedX = std::min(std::max(updatedPosition.x, 0.f), static_cast<float>(width - (2 * radius)));
    float boundedY = std::min(std::max(updatedPosition.y, 0.f), static_cast<float>(height - (2 * radius)));
    particle.setPosition(boundedX, boundedY);

    particleVelocity += gravity * elapsed.asSeconds();
}

float calculateDensity(const Vector2f &point) {
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

void visualizeDensity(sf::RenderWindow &window) {
    sf::Vector2u bounds = window.getSize();
    unsigned int width = bounds.x;
    unsigned int height = bounds.y;

    const float scale = 8.0;

    for (int i = 0; i < width / scale; i++) {
        for (int j = 0; j < height / scale; j++) {
            const Vector2f point(i * scale, j * scale);
            const float density = calculateDensity(point);
            sf::RectangleShape pixel(Vector2f(scale, scale));
            pixel.setPosition(Vector2f(i*scale, j*scale));
            pixel.setFillColor(sf::Color(0, 0, 255, 255 * std::min(density * 5e2, 1.0)));
            window.draw(pixel);
        }
    }
}

int main()
{
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;

    auto window = sf::RenderWindow{ { 1280u, 720u }, "Fluid Sim Project" , sf::Style::Default, settings};
    window.setFramerateLimit(60);

    // get the size of the window
    sf::Vector2u bounds = window.getSize();
    unsigned int width = bounds.x;
    unsigned int height = bounds.y;

    float radius = 2.0;

    // spawnParticles(50, 50, radius, Vector2f(static_cast<float>(width)/2.f, static_cast<float>(height)/2.f), 2);
    spawnRandomParticles(1000, radius, bounds);

    sf::Clock clock;
    Vector2f center(static_cast<float>(width)/2.f, static_cast<float>(height)/2.f);
    Circle densityCircle(center, smoothingRadius);
    densityCircle.setFillColor(sf::Color(0, 0, 0, 0));
    densityCircle.setOutlineColor(sf::Color(0, 255, 255));
    densityCircle.setOutlineThickness(1.0);

    float density = calculateDensity(center);
    std::cout << "density: " << density << std::endl;

    bool leftMousePressed = false;
    bool rightMousePressed = false;
    // float prevMouse[2];
    Vector2f prevMouse(0.f, 0.f);

    while (window.isOpen())
    {
        for (auto event = sf::Event{}; window.pollEvent(event);)
        {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left)
                {
                    leftMousePressed = true;
                    std::cout << "left mouse pressed" << std::endl;

                    center = Vector2f(event.mouseButton.x, event.mouseButton.y);
                    densityCircle.setPosition(center - Vector2f(smoothingRadius, smoothingRadius));
                    float density = calculateDensity(center);
                    std::cout << "density: " << density << std::endl;
                }
                if (event.mouseButton.button == sf::Mouse::Right)
                {
                    rightMousePressed = true;
                    prevMouse.x = event.mouseButton.x;
                    prevMouse.y = event.mouseButton.y;
                }
            }
            if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left)
                {
                    leftMousePressed = false;
                    std::cout << "left mouse released" << std::endl;

                    center = Vector2f(event.mouseButton.x, event.mouseButton.y);
                    densityCircle.setPosition(center - Vector2f(smoothingRadius, smoothingRadius));
                    float density = calculateDensity(center);
                    std::cout << "density: " << density << std::endl;
                }
                if (event.mouseButton.button == sf::Mouse::Right)
                {
                    rightMousePressed = false;
                }
            }

            if (event.type == sf::Event::MouseMoved)
            {
                if (leftMousePressed)
                {
                    center = Vector2f(event.mouseMove.x, event.mouseMove.y);
                    densityCircle.setPosition(center - Vector2f(smoothingRadius, smoothingRadius));
                }

                if (rightMousePressed)
                {
                    Vector2f centerDiff = Vector2f(prevMouse - center).normalized();

                    Vector2f mousePos(event.mouseMove.x, event.mouseMove.y);
                    Vector2f mouseDiff = mousePos - prevMouse;

                    float mouseToCenterProj = mouseDiff.inner_product(centerDiff);
                    smoothingRadius = std::max(smoothingRadius + mouseToCenterProj, minSmoothingRadius);
                    // NaN check
                    if (smoothingRadius != smoothingRadius) {
                        smoothingRadius = minSmoothingRadius;
                    }
                    std::cout << smoothingRadius << std::endl;
                    densityCircle.setRadius(smoothingRadius);
                    densityCircle.setPosition(center - Vector2f(smoothingRadius, smoothingRadius));
                }

                prevMouse.x = event.mouseMove.x;
                prevMouse.y = event.mouseMove.y;

                if (leftMousePressed || rightMousePressed)
                {
                    float density = calculateDensity(center);
                    std::cout << "density: " << density << std::endl;
                }
            }
        }

        window.clear();

        // visualizeDensity(window);

        
        window.draw(densityCircle);

        sf::Time elapsed = clock.restart();
        for (int i = 0; i < particles.size(); i++) {
            Circle &particle = particles[i];
            Vector2f &particleVelocity = particleVelocities[i];

            // applyGravity(particle, particleVelocity, elapsed, bounds);

            window.draw(particle);
        }
        
        window.display();
    }
}