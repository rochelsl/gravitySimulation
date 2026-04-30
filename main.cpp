#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>

const int width = 1920;
const int height = 1080;

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;

    float radius;
};

int main() {
    auto window = sf::RenderWindow{{width, height}, "Particle Simulation"}; 
    window.setFramerateLimit(144); 

    std::vector<Particle> particles;

    const int N = 100;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::uniform_real_distribution<float> vel(-200.f, 200.f);
    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    float dx = static_cast<float>(width) / cols;
    float dy = static_cast<float>(height) / rows;

    for (int i = 0; i < N; ++i) {
        int ix = i % cols;
        int iy = i / cols;

        Particle p;
        p.radius = 8.f;
        p.position = {(ix + 0.5f) * dx + jitter(rng), (iy + 0.5f) * dy + jitter(rng)};
        p.velocity = {vel(rng), vel(rng)};
        p.acceleration = {0.f, 0.f};
        particles.push_back(p);
    }

    sf::CircleShape shape;
    shape.setFillColor(sf::Color::White);

    const float dt = 0.001f;

    while (window.isOpen()) { 
        for (auto event = sf::Event{}; window.pollEvent(event);) { 
            if (event.type == sf::Event::Closed) { 
                window.close(); 
            }
			else if (event.type == sf::Event::KeyPressed) {
				if (event.key.code == sf::Keyboard::Escape) {
					window.close();
				}
			}
        }
        window.clear();
        for (const auto& p : particles) {
            shape.setRadius(p.radius);
            shape.setOrigin({p.radius, p.radius});
            shape.setPosition(p.position);
            window.draw(shape);
        }
        window.display();
    }

    return 0;
}