#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

const int width = 1500;
const int height = 1500;

const float gConst = 1000.0;
const float particleMass = 100.0;

//For Barnes-Hut algorithm, to stop subdivision at some spatial resolutions to prevent exploding tree depth
const float MIN_SIZE = 1.0f;

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

bool contains(const Quad& q, const sf::Vector2f& pos) {
    return (pos.x >= q.x - q.halfSize &&
            pos.x <= q.x + q.halfSize &&
            pos.y >= q.y - q.halfSize &&
            pos.y <= q.y + q.halfSize);
}

//Inserting particles into the quadrant
void insert(Node* node, Particle* p) {
    //Deactivated for periodic boundary
    //if (!contains(node->boundary, p->position)) return;

    // existing code continues
    if (node->particle == nullptr && node->isLeaf()) {
        node->particle = p;
        return;
    }
    // empty leaf
    if (node->particle == nullptr && node->isLeaf()) {
        node->particle = p;
        return;
    }

    // subdivide if needed
    if (node->isLeaf()) {
        // stop subdivision if too small
        if (node->boundary.halfSize < MIN_SIZE) {
            // fallback: accumulate mass instead of subdividing
            if (node->particle == nullptr) {
                node->particle = p;
            }
            return;
        }

        subdivide(node);

        Particle* old = node->particle;
        node->particle = nullptr;

        if (old) {
            int qOld = getQuadrant(node->boundary, old->position);
            insert(node->children[qOld], old);
        }
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

    // periodic wrapping
    if (r.x >  width * 0.5f) r.x -= width;
    if (r.x < -width * 0.5f) r.x += width;
    if (r.y >  height * 0.5f) r.y -= height;
    if (r.y < -height * 0.5f) r.y += height;

    const float eps = 5.0f;
    float dist2 = r.x * r.x + r.y * r.y + eps * eps;
    float dist = std::sqrt(dist2);

    float s = node->boundary.halfSize * 2.f;

    // ✅ ONLY here decide: approximate OR recurse
    if (node->isLeaf() || (s / dist) < theta) {
        if (node->particle == &p) return;

        float invDist = 1.0f / dist;
        float invDist3 = invDist * invDist * invDist;

        float f = gConst * node->mass * invDist3;
        p.acceleration += r * f;
    } else {
        for (int i = 0; i < 4; ++i) {
            if (node->children[i])
                computeForce(node->children[i], p, theta);
        }
    }

    // optional repulsion (keep small)
    const float rCut = 10.0f;
    if (dist < rCut) {
        float rep = 0.05f * (rCut - dist);
        p.acceleration -= (r / dist) * rep;
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
    root->boundary = {
    width / 2.f,
    height / 2.f,
    static_cast<float>(std::max(width, height))* 0.5f
};
    //root->boundary = {width / 2.f, height / 2.f, static_cast<float>(std::max(width, height))};

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

    //perdiodic boundary
    for (auto& p : particles) {
    if (p.position.x < 0) p.position.x += width;
    if (p.position.x >= width) p.position.x -= width;

    if (p.position.y < 0) p.position.y += height;
    if (p.position.y >= height) p.position.y -= height;
}
}

int main() {
    auto window = sf::RenderWindow{{width, height}, "Particle Simulation"}; 
    window.setFramerateLimit(144); 

    std::vector<Particle> particles;

    const int N = 5000;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::uniform_real_distribution<float> vel(-20.f, 20.f);
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

    //Remove net momentum which might occur from random initial velocity initialization
    sf::Vector2f vcm = {0.f, 0.f};
    for (const auto& p : particles)
        vcm += p.velocity;
    vcm /= static_cast<float>(particles.size());
    for (auto& p : particles)
        p.velocity -= vcm;

    // Initial acceleration for velocity Verlet.
    computeGravityBH(particles);

    sf::CircleShape shape;
    shape.setFillColor(sf::Color::White);

    const float dt = 0.0005f;

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
        for (auto& p : particles) {
            p.velocity += 0.5f * p.acceleration * dt;
            //optional cooling/dissipative term so that stable clusters are formed
            p.velocity *= 0.9999f;
        }

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