#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>

const int width = 1920;
const int height = 1080;

const float gConst = 1000;

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;

    float mass;
    float radius;
};

void computeGravity(std::vector<Particle>& particles) {
    // Reset accelerations 
    for (auto& p : particles) {
        p.acceleration = {0.f, 0.f};
    }

    for (Particle pi: particles) {
        for (Particle pj: particles) {
            sf::Vector2f r = pi.position - pj.position;
            float r2 = r.x * r.x + r.y * r.y;

            float f = gConst * (pi.mass * pj.mass) / r2;

            sf::Vector2f force = r * f;

            pi.acceleration -= force;
            pj.acceleration += force;

        }
    }
}

void integrate(std::vector<Particle>& particles, float dt) {
    for (auto& p : particles) {
        p.position += p.velocity * dt + 0.5f * p.acceleration * dt * dt;
        p.velocity += 0.5f * p.acceleration * dt;
    }
}

int main() {
    auto window = sf::RenderWindow{{width, height}, "Particle Simulation"}; 
    window.setFramerateLimit(144); 

    std::vector<Particle> particles;

    const int N = 1000;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::uniform_real_distribution<float> vel(-10.f, 10.f);
    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    float dx = static_cast<float>(width) / cols;
    float dy = static_cast<float>(height) / rows;

    for (int i = 0; i < N; ++i) {
        int ix = i % cols;
        int iy = i / cols;

        Particle p;
        p.radius = 3.f;
        p.position = {(ix + 0.5f) * dx + jitter(rng), (iy + 0.5f) * dy + jitter(rng)};
        p.velocity = {vel(rng), vel(rng)};
        p.acceleration = {0.f, 0.f};
        particles.push_back(p);
    }

    // Initial acceleration for velocity Verlet.
    computeGravity(particles);

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

        // 1. integrate positions + half velocity
        integrate(particles, dt);
        // 2. recompute forces
        computeGravity(particles);
        // 3. finish velocity update
        for (auto& p : particles)
            p.velocity += 0.5f * p.acceleration * dt;

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