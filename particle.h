#pragma once

#include <SFML/Graphics.hpp>

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;
    float mass;
    float radius;
};

struct GPUPosition {
    float x, y, z, mass;
};

struct GPUVelocity {
    float vx, vy, vz, radius;
};

struct GPUAcceleration {
    float ax, ay, az, unused;
};

struct GPUForce {
    float ax, ay, unused1, unused2;
};
