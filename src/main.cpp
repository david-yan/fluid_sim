#include <SFML/Graphics.hpp>
#include <iostream>

std::vector<sf::CircleShape> particles;
std::vector<sf::Vector2f> particleVelocities;

void spawnParticles(const int rows, const int cols, const float radius, const sf::Vector2f &center) {
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
            sf::CircleShape particle(radius);
            particle.setFillColor(sf::Color(0, 0, 255));

            sf::Vector2f displacement((i - rowCenterIdx) * 2 * radius, (j - colCenterIdx) * 2 * radius);
            particle.setPosition(center + displacement);
            particles.push_back(particle);
            particleVelocities.push_back(sf::Vector2f(0.f, 0.f));
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
    sf::Vector2u size = window.getSize();
    unsigned int width = size.x;
    unsigned int height = size.y;

    float circleRadius = 5;
    // sf::CircleShape circle(circleRadius);
    // // set the shape color to blue
    // circle.setFillColor(sf::Color(0, 0, 255));
    // // set the number of sides (quality)
    // // circle.setPointCount(100);
    // circle.setPosition(width/2, height/2);
    // sf::Vector2f circleVelocity(5.f, 0.f);

    spawnParticles(10, 10, circleRadius, sf::Vector2f(float(width)/2.f, float(height)/2.f));

    sf::Vector2f gravity(0.f, 9.8);
    sf::Vector2f damping(0.9, 0.9);

    sf::Clock clock;

    while (window.isOpen())
    {
        for (auto event = sf::Event{}; window.pollEvent(event);)
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
        }

        window.clear();

        sf::Time elapsed = clock.restart();
        for (int i = 0; i < particles.size(); i++) {
            sf::CircleShape &particle = particles[i];
            sf::Vector2f &particleVelocity = particleVelocities[i];

            particle.move(particleVelocity);

            // Check particle bounds
            sf::Vector2f updatedPosition = particle.getPosition();
            if (updatedPosition.x < 0 || updatedPosition.x + (2 * circleRadius) >= width) {
                particleVelocity.x *= -1 * damping.x;
            }
            if (updatedPosition.y < 0 || updatedPosition.y + (2 * circleRadius) >= height) {
                particleVelocity.y *= -1 * damping.y;
            }
            float boundedX = std::min(std::max(updatedPosition.x, 0.f), float(width - (2 * circleRadius)));
            float boundedY = std::min(std::max(updatedPosition.y, 0.f), float(height - (2 * circleRadius)));
            particle.setPosition(boundedX, boundedY);

            particleVelocity += gravity * elapsed.asSeconds();

            window.draw(particle);
        }
        
        window.display();
    }
}