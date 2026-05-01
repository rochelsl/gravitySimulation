//When eomploying OpenGL run using: g++ -std=c++17 main.cpp -o main -lGLEW -lGL -lGLU -lsfml-graphics -lsfml-window -lsfml-system
#include <GL/glew.h>
#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <string>

const int width = 1500;
const int height = 1500;

const float gConst = 1000.0;
const float particleMass = 10.0;

//For Barnes-Hut algorithm, to stop subdivision at some spatial resolutions to prevent exploding tree depth
const float MIN_SIZE = 1.0f;

//Adding grid buffers for GPU computation
//Higher grid resolution e.g. 256, 144 gices sharper gravity but costs more
//Lower gives smoother, faster behavior
const int gridW = 32;
const int gridH = 18;
const int gridSize = gridW * gridH;

struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;

    float mass;
    float radius;
};

//Replacing Particle as the active simulation state
struct GPUPosition {
    float x, y, z, mass;
};

struct GPUVelocity {
    float vx, vy, vz, radius;
};

//Adding GPU acceleartion structs to prepare the Barnes-Hut algorithm for GPU computation
struct GPUAcceleration {
    float ax, ay, az, unused;
};

struct GPUForce {
    float ax, ay, unused1, unused2;
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

    // ONLY here decide: approximate OR recurse
    if (node->isLeaf()) {
        if (node->particle == nullptr || node->particle == &p) return;

        // optional repulsion (keep small)
        const float rCut = 10.0f;
        if (dist < rCut) {
            float rep = 0.05f * (rCut - dist);
            p.acceleration -= (r / dist) * rep;
        }

        float invDist = 1.0f / dist;
        float invDist3 = invDist * invDist * invDist;
        float f = gConst * node->mass * invDist3;
        p.acceleration += r * f;
    }
    else if ((s / dist) < theta) {
        float invDist = 1.0f / dist;
        float invDist3 = invDist * invDist * invDist;
        float f = gConst * node->mass * invDist3;
        p.acceleration += r * f;
    }
    else {
        for (int i = 0; i < 4; ++i)
            if (node->children[i])
                computeForce(node->children[i], p, theta);
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

//for virialized initial conditions, where 2T = |U| (T: total kinetic energy, U: total potential energy)
float computeKinetic(const std::vector<Particle>& particles) {
    float T = 0.f;
    for (const auto& p : particles) {
        float v2 = p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y;
        T += 0.5f * p.mass * v2;
    }
    return T;
}

float computePotential(const std::vector<Particle>& particles) {
    float U = 0.f;

    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            sf::Vector2f r = particles[j].position - particles[i].position;

            // periodic wrapping
            if (r.x >  width * 0.5f) r.x -= width;
            if (r.x < -width * 0.5f) r.x += width;
            if (r.y >  height * 0.5f) r.y -= height;
            if (r.y < -height * 0.5f) r.y += height;

            float dist = std::sqrt(r.x*r.x + r.y*r.y + 25.0f); // eps^2 = 5^2
            U -= gConst * particles[i].mass * particles[j].mass / dist;
        }
    }

    return U;
}

//GPU shader helper functions
GLuint compileShader(GLenum type, const char* source) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    GLint success = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    if (!success) {
        char log[1024];
        glGetShaderInfoLog(shader, 1024, nullptr, log);
        std::cerr << "Shader compile error:\n" << log << std::endl;
    }

    return shader;
}

GLuint createProgram(const char* vertexSrc, const char* fragmentSrc) {
    GLuint vs = compileShader(GL_VERTEX_SHADER, vertexSrc);
    GLuint fs = compileShader(GL_FRAGMENT_SHADER, fragmentSrc);

    GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    GLint success = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &success);

    if (!success) {
        char log[1024];
        glGetProgramInfoLog(program, 1024, nullptr, log);
        std::cerr << "Program link error:\n" << log << std::endl;
    }

    glDeleteShader(vs);
    glDeleteShader(fs);

    return program;
}

//Adding a compute program helper
GLuint createComputeProgram(const char* computeSrc) {
    GLuint cs = compileShader(GL_COMPUTE_SHADER, computeSrc);

    GLuint program = glCreateProgram();
    glAttachShader(program, cs);
    glLinkProgram(program);

    GLint success = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &success);

    if (!success) {
        char log[1024];
        glGetProgramInfoLog(program, 1024, nullptr, log);
        std::cerr << "Compute program link error:\n" << log << std::endl;
    }

    glDeleteShader(cs);
    return program;
}

int main() {
    sf::ContextSettings settings;
    settings.majorVersion = 4;
    settings.minorVersion = 3;
    settings.depthBits = 24;
    settings.stencilBits = 8;

    auto window = sf::RenderWindow(
        sf::VideoMode(width, height),
        "Particle Simulation",
        sf::Style::Default,
        settings
    );

    window.setActive(true);
    //Initialize GLEW after window creation
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";
        return -1;
    }

    const char* vertexShaderSrc = R"(
    #version 430 core
    layout(std430, binding = 0) buffer Positions {
        vec4 posMass[];
    };
    uniform float uWidth;
    uniform float uHeight;
    uniform float uPointSize;
    void main() {
        vec2 pos = posMass[gl_VertexID].xy;
        float x = (pos.x / uWidth) * 2.0 - 1.0;
        float y = 1.0 - (pos.y / uHeight) * 2.0;
        gl_Position = vec4(x, y, 0.0, 1.0);
        gl_PointSize = uPointSize;
    }
    )";
    const char* fragmentShaderSrc = R"(
    #version 430 core
    out vec4 FragColor;
    void main() {
        FragColor = vec4(1.0, 1.0, 1.0, 1.0);
    }
    )";

    GLuint renderProgram = createProgram(vertexShaderSrc, fragmentShaderSrc);

    //Adding the OpenGL shaders after GLEW initialization
    const char* integrateShaderSrc = R"(
    #version 430 core
    layout(local_size_x = 256) in;
    layout(std430, binding = 0) buffer Positions {
        vec4 posMass[];
    };
    layout(std430, binding = 1) buffer Velocities {
        vec4 velRadius[];
    };
    layout(std430, binding = 2) buffer Accelerations {
        vec4 acc[];
    };
    uniform float dt;
    uniform float width;
    uniform float height;
    uniform float damping;
    uniform uint numParticles;
    void main() {
        uint i = gl_GlobalInvocationID.x;
        if (i >= numParticles) return;
            vec2 pos = posMass[i].xy;
            vec2 vel = velRadius[i].xy;
            vec2 a   = acc[i].xy;

            vel += a * dt;
            vel *= damping;
            pos += vel * dt;
        // periodic boundary
        if (pos.x < 0.0) pos.x += width;
        if (pos.x >= width) pos.x -= width;
        if (pos.y < 0.0) pos.y += height;
        if (pos.y >= height) pos.y -= height;
        posMass[i].xy = pos;
        velRadius[i].xy = vel;
    }
    )";
    GLuint integrateProgram = createComputeProgram(integrateShaderSrc);

    //After creating renderProgram
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    //After integration shader
    const char* gravityShaderSrc = R"(
    #version 430 core
    layout(local_size_x = 256) in;
    layout(std430, binding = 0) buffer Positions {
        vec4 posMass[];
    };
    layout(std430, binding = 2) buffer Accelerations {
        vec4 acc[];
    };
    uniform float width;
    uniform float height;
    uniform float gConst;
    uniform float eps;
    uniform uint numParticles;
    void main() {
        uint i = gl_GlobalInvocationID.x;
        if (i >= numParticles) return;
        vec2 posI = posMass[i].xy;
        float massI = posMass[i].w;
        vec2 a = vec2(0.0);
        for (uint j = 0; j < numParticles; ++j) {
            if (i == j) continue;
            vec2 r = posMass[j].xy - posI;
            // periodic minimum-image convention
            if (r.x >  width * 0.5) r.x -= width;
            if (r.x < -width * 0.5) r.x += width;
            if (r.y >  height * 0.5) r.y -= height;
            if (r.y < -height * 0.5) r.y += height;
            float dist2 = dot(r, r) + eps * eps;
            float invDist = inversesqrt(dist2);
            float invDist3 = invDist * invDist * invDist;
            // acceleration, not force
            a += gConst * posMass[j].w * r * invDist3;
        }
        acc[i] = vec4(a, 0.0, 0.0);
    }
    )";
    GLuint gravityProgram = createComputeProgram(gravityShaderSrc);

    const char* clearGridShaderSrc = R"(
    #version 430 core
    layout(local_size_x = 256) in;
    layout(std430, binding = 3) buffer DensityGrid {
        uint density[];
    };
    uniform uint gridSize;
    void main() {
        uint i = gl_GlobalInvocationID.x;
        if (i >= gridSize) return;
        density[i] = 0u;
    }
    )";
    GLuint clearGridProgram = createComputeProgram(clearGridShaderSrc);

    const char* depositShaderSrc = R"(
    #version 430 core
    layout(local_size_x = 256) in;
    layout(std430, binding = 0) buffer Positions {
        vec4 posMass[];
    };
    layout(std430, binding = 3) buffer DensityGrid {
        uint density[];
    };
    uniform float width;
    uniform float height;
    uniform uint gridW;
    uniform uint gridH;
    uniform uint numParticles;
    void main() {
        uint i = gl_GlobalInvocationID.x;
        if (i >= numParticles) return;
        vec2 pos = posMass[i].xy;
        uint cellX = uint(clamp(pos.x / width  * float(gridW), 0.0, float(gridW - 1)));
        uint cellY = uint(clamp(pos.y / height * float(gridH), 0.0, float(gridH - 1)));
        uint cellIndex = cellY * gridW + cellX;
        atomicAdd(density[cellIndex], 1u);
    }
    )";
    GLuint depositProgram = createComputeProgram(depositShaderSrc);

    const char* computeGridForceShaderSrc = R"(
    #version 430 core
    layout(local_size_x = 256) in;
    layout(std430, binding = 3) buffer DensityGrid {
        uint density[];
    };
    layout(std430, binding = 4) buffer ForceGrid {
        vec4 force[];
    };
    uniform uint gridW;
    uniform uint gridH;
    uniform float width;
    uniform float height;
    uniform float gConst;
    uniform float particleMass;
    uniform float eps;
    void main() {
        uint idx = gl_GlobalInvocationID.x;
        uint gridSize = gridW * gridH;
        if (idx >= gridSize) return;
        uint x = idx % gridW;
        uint y = idx / gridW;
        vec2 cellSize = vec2(width / float(gridW), height / float(gridH));
        vec2 posA = vec2(
            (float(x) + 0.5) * cellSize.x,
            (float(y) + 0.5) * cellSize.y
        );
        vec2 acc = vec2(0.0);
        for (uint j = 0; j < gridSize; ++j) {
            uint count = density[j];
            if (count == 0u || j == idx) continue;
            uint xj = j % gridW;
            uint yj = j / gridW;
            vec2 posB = vec2(
                (float(xj) + 0.5) * cellSize.x,
                (float(yj) + 0.5) * cellSize.y
            );
            vec2 r = posB - posA;
            // periodic minimum-image convention
            if (r.x >  width * 0.5) r.x -= width;
            if (r.x < -width * 0.5) r.x += width;
            if (r.y >  height * 0.5) r.y -= height;
            if (r.y < -height * 0.5) r.y += height;
            float mass = float(count) * particleMass;
            float dist2 = dot(r, r) + eps * eps;
            float invDist = inversesqrt(dist2);
            float invDist3 = invDist * invDist * invDist;
            acc += gConst * mass * r * invDist3;
        }
        force[idx] = vec4(acc, 0.0, 0.0);
    }
    )";
    GLuint computeGridForceProgram = createComputeProgram(computeGridForceShaderSrc);

    const char* sampleForceShaderSrc = R"(
    #version 430 core
    layout(local_size_x = 256) in;
    layout(std430, binding = 0) buffer Positions {
        vec4 posMass[];
    };
    layout(std430, binding = 2) buffer Accelerations {
        vec4 acc[];
    };
    layout(std430, binding = 4) buffer ForceGrid {
        vec4 force[];
    };
    uniform float width;
    uniform float height;
    uniform uint gridW;
    uniform uint gridH;
    uniform uint numParticles;
    void main() {
        uint i = gl_GlobalInvocationID.x;
        if (i >= numParticles) return;
        vec2 pos = posMass[i].xy;
        uint cellX = uint(clamp(pos.x / width  * float(gridW), 0.0, float(gridW - 1)));
        uint cellY = uint(clamp(pos.y / height * float(gridH), 0.0, float(gridH - 1)));
        uint cellIndex = cellY * gridW + cellX;
        acc[i] = force[cellIndex];
    }
    )";
    GLuint sampleForceProgram = createComputeProgram(sampleForceShaderSrc);

    std::vector<Particle> particles;

    const int N = 30000;
    particles.reserve(N);

    std::mt19937 rng(std::random_device{}());

    std::uniform_real_distribution<float> vel(-100.f, 100.f);
    // Grid initialization avoids catastrophic LJ overlaps from random placement.
    int cols = static_cast<int>(std::ceil(std::sqrt(N * static_cast<float>(width) / height)));
    int rows = static_cast<int>(std::ceil(static_cast<float>(N) / cols));
    float dx = static_cast<float>(width) / cols;
    float dy = static_cast<float>(height) / rows;
    //deviations from perfect square lattice
    std::uniform_real_distribution<float> jitterX(-0.8f * dx, 0.8f * dx);
    std::uniform_real_distribution<float> jitterY(-0.8f * dy, 0.8f * dy);

    for (int i = 0; i < N; ++i) {
        int ix = i % cols;
        int iy = i / cols;

        Particle p;
        p.radius = 1.f;
        p.position = {
            (ix + 0.5f) * dx + jitterX(rng),
            (iy + 0.5f) * dy + jitterY(rng)
        };
        p.velocity = {vel(rng), vel(rng)};
        p.acceleration = {0.f, 0.f};
        p.mass = particleMass;
        particles.push_back(p);
    }

    //virialized initial conditions to rescale velocities to virial equilibrium
    float T = computeKinetic(particles);
    float U = computePotential(particles);

    // target: 2T = |U|
    //if the system is too "hot"/gas-like, reduce the scaling factor
    float scale = 0.01f * std::sqrt(std::abs(U) / (2.0f * T));

    // rescale velocities
    for (auto& p : particles) {
        p.velocity *= scale;
    }

    //Add a rotational velocity around the box center
    sf::Vector2f center = {
    width * 0.5f,
    height * 0.5f
    };
    const float rotationStrength = 100.0f; // tune this
    for (auto& p : particles) {
        sf::Vector2f r = p.position - center;
        // periodic minimum-image displacement from center
        if (r.x >  width * 0.5f) r.x -= width;
        if (r.x < -width * 0.5f) r.x += width;
        if (r.y >  height * 0.5f) r.y -= height;
        if (r.y < -height * 0.5f) r.y += height;
        float dist = std::sqrt(r.x * r.x + r.y * r.y) + 1e-6f;
        // tangential direction
        sf::Vector2f tangent = {
            -r.y / dist,
            r.x / dist
        };
        // optional: weaker rotation near center, stronger farther out
        float speed = rotationStrength * std::sqrt(dist / std::max(width, height));
        p.velocity += tangent * speed;
    }

    //Remove net momentum/drift which might occur from random initial velocity initialization
    sf::Vector2f vcm = {0.f, 0.f};
    for (const auto& p : particles)
        vcm += p.velocity;
    vcm /= static_cast<float>(particles.size());
    for (auto& p : particles)
        p.velocity -= vcm;

    //after initialization, virial scaling, rotation and drift removal add the GPU arrays
    std::vector<GPUPosition> gpuPositions(particles.size());
    std::vector<GPUVelocity> gpuVelocities(particles.size());
    GLuint posSSBO;
    GLuint velSSBO;
    for (size_t i = 0; i < particles.size(); ++i) {
        gpuPositions[i] = {
            particles[i].position.x,
            particles[i].position.y,
            0.0f,
            particles[i].mass
        };
    }

    glGenBuffers(1, &posSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, posSSBO);
    glBufferData(
        GL_SHADER_STORAGE_BUFFER,
        gpuPositions.size() * sizeof(GPUPosition),
        gpuPositions.data(),
        GL_DYNAMIC_DRAW
    );
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, posSSBO);

    for (size_t i = 0; i < particles.size(); ++i) {
        gpuVelocities[i] = {
            particles[i].velocity.x,
            particles[i].velocity.y,
            0.0f,
            particles[i].radius
        };
    }

    glGenBuffers(1, &velSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, velSSBO);
    glBufferData(
        GL_SHADER_STORAGE_BUFFER,
        gpuVelocities.size() * sizeof(GPUVelocity),
        gpuVelocities.data(),
        GL_DYNAMIC_DRAW
    );
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, velSSBO);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    //After creating gpuPositions and gpuVelocities
    std::vector<GPUAcceleration> gpuAccelerations(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
        gpuAccelerations[i] = {
            particles[i].acceleration.x,
            particles[i].acceleration.y,
            0.0f,
            0.0f
        };
    }

    //Create SSBO
    GLuint accSSBO;
    glGenBuffers(1, &accSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, accSSBO);
    glBufferData(
        GL_SHADER_STORAGE_BUFFER,
        gpuAccelerations.size() * sizeof(GPUAcceleration),
        gpuAccelerations.data(),
        GL_DYNAMIC_DRAW
    );
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, accSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    std::vector<unsigned int> densityGrid(gridSize, 0);
    GLuint densitySSBO;
    glGenBuffers(1, &densitySSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, densitySSBO);
    glBufferData(
        GL_SHADER_STORAGE_BUFFER,
        densityGrid.size() * sizeof(unsigned int),
        densityGrid.data(),
        GL_DYNAMIC_DRAW
    );
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, densitySSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    //After densitySSBO creation
    std::vector<GPUForce> forceGrid(gridSize);
    GLuint forceGridSSBO;
    glGenBuffers(1, &forceGridSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, forceGridSSBO);
    glBufferData(
        GL_SHADER_STORAGE_BUFFER,
        forceGrid.size() * sizeof(GPUForce),
        forceGrid.data(),
        GL_DYNAMIC_DRAW
    );
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, forceGridSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    const float dt = 0.0005f;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
			else if (event.type == sf::Event::KeyPressed) {
				if (event.key.code == sf::Keyboard::Escape) {
					window.close();
				}
			}
        }

        GLuint particleGroups = static_cast<GLuint>((particles.size() + 255) / 256);
        GLuint gridGroups = static_cast<GLuint>((gridSize + 255) / 256);

        // 1. Clear density grid
        glUseProgram(clearGridProgram);
        glUniform1ui(
            glGetUniformLocation(clearGridProgram, "gridSize"),
            static_cast<unsigned int>(gridSize)
        );
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, densitySSBO);
        glDispatchCompute(gridGroups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        // 2. Deposit particles into density grid
        glUseProgram(depositProgram);
        glUniform1f(glGetUniformLocation(depositProgram, "width"), static_cast<float>(width));
        glUniform1f(glGetUniformLocation(depositProgram, "height"), static_cast<float>(height));
        glUniform1ui(glGetUniformLocation(depositProgram, "gridW"), static_cast<unsigned int>(gridW));
        glUniform1ui(glGetUniformLocation(depositProgram, "gridH"), static_cast<unsigned int>(gridH));
        glUniform1ui(glGetUniformLocation(depositProgram, "numParticles"), static_cast<unsigned int>(particles.size()));
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, posSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, densitySSBO);
        glDispatchCompute(particleGroups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        // 3. Compute grid force field
        glUseProgram(computeGridForceProgram);
        glUniform1ui(glGetUniformLocation(computeGridForceProgram, "gridW"), static_cast<unsigned int>(gridW));
        glUniform1ui(glGetUniformLocation(computeGridForceProgram, "gridH"), static_cast<unsigned int>(gridH));
        glUniform1f(glGetUniformLocation(computeGridForceProgram, "width"), static_cast<float>(width));
        glUniform1f(glGetUniformLocation(computeGridForceProgram, "height"), static_cast<float>(height));
        glUniform1f(glGetUniformLocation(computeGridForceProgram, "gConst"), gConst);
        glUniform1f(glGetUniformLocation(computeGridForceProgram, "particleMass"), particleMass);
        glUniform1f(glGetUniformLocation(computeGridForceProgram, "eps"), 20.0f);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, densitySSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, forceGridSSBO);
        glDispatchCompute(gridGroups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        // 4. Sample grid force onto particles
        glUseProgram(sampleForceProgram);
        glUniform1f(glGetUniformLocation(sampleForceProgram, "width"), static_cast<float>(width));
        glUniform1f(glGetUniformLocation(sampleForceProgram, "height"), static_cast<float>(height));
        glUniform1ui(glGetUniformLocation(sampleForceProgram, "gridW"), static_cast<unsigned int>(gridW));
        glUniform1ui(glGetUniformLocation(sampleForceProgram, "gridH"), static_cast<unsigned int>(gridH));
        glUniform1ui(glGetUniformLocation(sampleForceProgram, "numParticles"), static_cast<unsigned int>(particles.size()));
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, posSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, accSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, forceGridSSBO);
        glDispatchCompute(particleGroups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        // 5. GPU integration
        glUseProgram(integrateProgram);
        glUniform1f(glGetUniformLocation(integrateProgram, "dt"), dt);
        glUniform1f(glGetUniformLocation(integrateProgram, "width"), static_cast<float>(width));
        glUniform1f(glGetUniformLocation(integrateProgram, "height"), static_cast<float>(height));
        glUniform1f(glGetUniformLocation(integrateProgram, "damping"), 0.99999f);
        glUniform1ui(glGetUniformLocation(integrateProgram, "numParticles"), static_cast<unsigned int>(particles.size()));
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, posSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, velSSBO);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, accSSBO);
        glDispatchCompute(particleGroups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT);

        glViewport(0, 0, width, height);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(renderProgram);
        glUniform1f(glGetUniformLocation(renderProgram, "uWidth"), static_cast<float>(width));
        glUniform1f(glGetUniformLocation(renderProgram, "uHeight"), static_cast<float>(height));
        glUniform1f(glGetUniformLocation(renderProgram, "uPointSize"), 1.5f);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, posSSBO);
        glBindVertexArray(vao);
        glDrawArrays(GL_POINTS, 0, particles.size());

        // static int frame = 0;
        // frame++;
        // if (frame % 300 == 0) {
        //     glBindBuffer(GL_SHADER_STORAGE_BUFFER, densitySSBO);
        //     glGetBufferSubData(
        //         GL_SHADER_STORAGE_BUFFER,
        //         0,
        //         densityGrid.size() * sizeof(unsigned int),
        //         densityGrid.data()
        //     );
        //     unsigned int total = 0;
        //     for (auto v : densityGrid)
        //         total += v;
        //     std::cout << "Deposited particles: " << total << std::endl;
        // }

        window.display();
    }
    return 0;
}