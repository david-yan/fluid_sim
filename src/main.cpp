#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>

#include "constants.hpp"
#include "circle.hpp"
#include "vector2f.hpp"
#include "timer.hpp"
#include "particles.hpp"

float minSmoothingRadius(1.0);
float smoothingRadius(50.0);


void visualizeDensity(sf::RenderWindow &window, ParticleStore &particles) {
    sf::Vector2u bounds = window.getSize();
    unsigned int width = bounds.x;
    unsigned int height = bounds.y;

    for (int i = 0; i < particles.allDensities.size(); i++) {
        for (int j = 0; j < particles.allDensities[0].size(); j++) {
            const float density = particles.allDensities[i][j];
            sf::RectangleShape pixel(Vector2f(DENSITY_SCALE, DENSITY_SCALE));
            pixel.setPosition(Vector2f(i*DENSITY_SCALE, j*DENSITY_SCALE));
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

    // initialize particles
    ParticleStore particles(bounds, smoothingRadius);
    float radius = 2.0;

    // particles.spawnParticles(50, 50, radius, Vector2f(static_cast<float>(width)/2.f, static_cast<float>(height)/2.f), 2);
    particles.spawnRandomParticles(1000, radius, bounds);

    sf::Clock clock;
    Vector2f center(static_cast<float>(width)/2.f, static_cast<float>(height)/2.f);
    Circle densityCircle(center, smoothingRadius);
    densityCircle.setFillColor(sf::Color(0, 0, 0, 0));
    densityCircle.setOutlineColor(sf::Color(0, 255, 255));
    densityCircle.setOutlineThickness(1.0);

    float density = particles.calculateDensity(center);
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
                    float density = particles.calculateDensity(center);
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
                    float density = particles.calculateDensity(center);
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
                    float density = particles.calculateDensity(center);
                    std::cout << "density: " << density << std::endl;
                }
            }
        }

        window.clear();

        visualizeDensity(window, particles);

        window.draw(densityCircle);

        // sf::Time elapsed = clock.restart();
        // particles.applyGravity(elapsed);

        particles.calculateAllDensities();
        particles.visualizeQuadrants(window);
        particles.drawParticles(window);
        
        window.display();
    }
}