#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

const int width = 1920;
const int height = 1080;

const float gConst = 10;
const float particleMass = 1000;

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;

    float mass;
    float radius;
};

//For implementation of a Barnes-Hut algorithm
struct Quad {
    float x, y;     // center
    float halfSize; // half-width
};

struct Node {
    Quad boundary;

    float mass = 0.f;
    sf::Vector2f com = {0.f, 0.f}; // center of mass

    Particle* particle = nullptr; // leaf node (1 particle max)

    Node* children[4] = {nullptr, nullptr, nullptr, nullptr};

    bool isLeaf() const {
        return children[0] == nullptr;
    }
};

//Quadrant indexing
int getQuadrant(const Quad& q, const sf::Vector2f& pos) {
    int index = 0;
    if (pos.x > q.x) index += 1;
    if (pos.y > q.y) index += 2;
    return index;
}

//Subdivision of quadrants
void subdivide(Node* node) {
    float hs = node->boundary.halfSize * 0.5f;

    for (int i = 0; i < 4; ++i) {
        node->children[i] = new Node();
        node->children[i]->boundary.halfSize = hs;

        node->children[i]->boundary.x = node->boundary.x + hs * (i % 2 ? 1.f : -1.f);
        node->children[i]->boundary.y = node->boundary.y + hs * (i / 2 ? 1.f : -1.f);
    }
}

//Inserting particles into the quadrant
void insert(Node* node, Particle* p) {
    // empty leaf
    if (node->particle == nullptr && node->isLeaf()) {
        node->particle = p;
        return;
    }

    // subdivide if needed
    if (node->isLeaf()) {
        subdivide(node);

        // reinsert existing particle
        Particle* old = node->particle;
        node->particle = nullptr;

        int qOld = getQuadrant(node->boundary, old->position);
        insert(node->children[qOld], old);
    }

    int q = getQuadrant(node->boundary, p->position);
    insert(node->children[q], p);
}

//Computing the mass distribution from the bottom-up
void computeMass(Node* node) {
    if (node->isLeaf()) {
        if (node->particle) {
            node->mass = node->particle->mass;
            node->com = node->particle->position;
        }
        return;
    }

    node->mass = 0.f;
    node->com = {0.f, 0.f};

    for (int i = 0; i < 4; ++i) {
        Node* child = node->children[i];
        if (!child) continue;

        computeMass(child);

        node->mass += child->mass;
        node->com += child->com * child->mass;
    }

    if (node->mass > 0.f)
        node->com /= node->mass;
}

//Barnes-Hut force computation
void computeForce(Node* node, Particle& p, float theta) {
    if (node->mass == 0.f) return;

    sf::Vector2f r = node->com - p.position;
    float dist2 = r.x * r.x + r.y * r.y + 1e-6f;
    float dist = std::sqrt(dist2);

    // region size
    float s = node->boundary.halfSize * 2.f;

    // if far enough → approximate
    if (node->isLeaf() || (s / dist) < theta) {
        if (node->particle == &p) return; // skip self

        float f = gConst * (p.mass * node->mass) / dist2;
        sf::Vector2f force = (r / dist) * f;

        p.acceleration += force / p.mass;
    } else {
        // recurse
        for (int i = 0; i < 4; ++i) {
            if (node->children[i])
                computeForce(node->children[i], p, theta);
        }
    }
}

//For memory managment, the tree each frame has to be deleted
void deleteTree(Node* node) {
    if (!node) return;

    for (int i = 0; i < 4; ++i)
        deleteTree(node->children[i]);

    delete node;
}

void computeGravityBH(std::vector<Particle>& particles) {
    // reset accelerations
    for (auto& p : particles)
        p.acceleration = {0.f, 0.f};

    // build root node (cover entire domain)
    Node* root = new Node();
    root->boundary = {width / 2.f, height / 2.f, static_cast<float>(std::max(width, height))};

    // insert particles
    for (auto& p : particles)
        insert(root, &p);

    // compute mass distribution
    computeMass(root);

    // compute forces
    float theta = 0.5f;

    for (auto& p : particles)
        computeForce(root, p, theta);

    //Freeing memory after each force computation step/frame
    deleteTree(root);
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

    const int N = 10000;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::uniform_real_distribution<float> vel(-0.f, 0.f);
    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    float dx = static_cast<float>(width) / cols;
    float dy = static_cast<float>(height) / rows;

    for (int i = 0; i < N; ++i) {
        int ix = i % cols;
        int iy = i / cols;

        Particle p;
        p.radius = 1.f;
        p.position = {(ix + 0.5f) * dx + jitter(rng), (iy + 0.5f) * dy + jitter(rng)};
        p.velocity = {vel(rng), vel(rng)};
        p.acceleration = {0.f, 0.f};
        p.mass = particleMass;
        particles.push_back(p);
    }

    // Initial acceleration for velocity Verlet.
    computeGravityBH(particles);

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
        computeGravityBH(particles);
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

        // float ts = particles[0].position.x;
        // std::cout << ts << std::endl;
    }


    return 0;
}