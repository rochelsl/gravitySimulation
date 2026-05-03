#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <string>

#include "constants.h"

namespace {

float computeKinetic(const std::vector<Particle>& particles) {
    float t = 0.f;
    for (const auto& p : particles) {
        float v2 = p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y;
        t += 0.5f * p.mass * v2;
    }
    return t;
}

float computePotential(const std::vector<Particle>& particles) {
    float u = 0.f;
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            sf::Vector2f r = particles[j].position - particles[i].position;
            if (r.x > kWidth * 0.5f) r.x -= kWidth;
            if (r.x < -kWidth * 0.5f) r.x += kWidth;
            if (r.y > kHeight * 0.5f) r.y -= kHeight;
            if (r.y < -kHeight * 0.5f) r.y += kHeight;
            float dist = std::sqrt(r.x * r.x + r.y * r.y + 25.0f);
            u -= kGConst * particles[i].mass * particles[j].mass / dist;
        }
    }
    return u;
}

}

std::vector<Particle> createParticles(int count) {
    std::vector<Particle> particles;
    particles.reserve(count);

    std::mt19937 rng(std::random_device{}());

    if (std::string(kInitType) == "grid") {
        int cols = static_cast<int>(std::ceil(std::sqrt(count * static_cast<float>(kWidth) / kHeight)));
        int rows = static_cast<int>(std::ceil(static_cast<float>(count) / cols));
        float dx = static_cast<float>(kWidth) / cols;
        float dy = static_cast<float>(kHeight) / rows;
        std::uniform_real_distribution<float> jitterX(-0.8f * dx, 0.8f * dx);
        std::uniform_real_distribution<float> jitterY(-0.8f * dy, 0.8f * dy);
        std::uniform_real_distribution<float> vel(-100.f, 100.f);
        for (int i = 0; i < count; ++i) {
            int ix = i % cols;
            int iy = i / cols;
            Particle p{};
            p.radius = 1.f;
            p.position = {(ix + 0.5f) * dx + jitterX(rng), (iy + 0.5f) * dy + jitterY(rng)};
            p.velocity = {vel(rng), vel(rng)};
            p.acceleration = {0.f, 0.f};
            p.mass = kParticleMass;
            particles.push_back(p);
        }
    } else {
        float widthNorm1 = kWidth * 0.5f;
        float widthNorm2 = kWidth * 0.25f;
        float heightNorm1 = kHeight * 0.5f;
        float heightNorm2 = kHeight * 0.25f;
        std::normal_distribution<float> posX(widthNorm1, widthNorm2);
        std::normal_distribution<float> posY(heightNorm1, heightNorm2);
        std::uniform_real_distribution<float> vel(-5.f, 5.f);

        for (int i = 0; i < count; ++i) {
            Particle p{};
            p.radius = 1.f;
            float x = posX(rng);
            float y = posY(rng);
            while (x < 0) x += kWidth;
            while (x >= kWidth) x -= kWidth;
            while (y < 0) y += kHeight;
            while (y >= kHeight) y -= kHeight;
            p.position = {x, y};
            p.velocity = {vel(rng), vel(rng)};
            p.acceleration = {0.f, 0.f};
            p.mass = kParticleMass;
            particles.push_back(p);
        }
    }

    float t = computeKinetic(particles);
    float u = computePotential(particles);
    float scale = 0.01f * std::sqrt(std::abs(u) / (2.0f * t));
    for (auto& p : particles) p.velocity *= scale;

    sf::Vector2f center = {kWidth * 0.5f, kHeight * 0.5f};
    constexpr float rotationStrength = 10.0f;
    for (auto& p : particles) {
        sf::Vector2f r = p.position - center;
        if (r.x > kWidth * 0.5f) r.x -= kWidth;
        if (r.x < -kWidth * 0.5f) r.x += kWidth;
        if (r.y > kHeight * 0.5f) r.y -= kHeight;
        if (r.y < -kHeight * 0.5f) r.y += kHeight;
        float dist = std::sqrt(r.x * r.x + r.y * r.y) + 1e-6f;
        sf::Vector2f tangent = {-r.y / dist, r.x / dist};
        float speed = rotationStrength * std::sqrt(dist / std::max(kWidth, kHeight));
        p.velocity += tangent * speed;
    }

    sf::Vector2f vcm = {0.f, 0.f};
    for (const auto& p : particles) vcm += p.velocity;
    vcm /= static_cast<float>(particles.size());
    for (auto& p : particles) p.velocity -= vcm;

    return particles;
}
