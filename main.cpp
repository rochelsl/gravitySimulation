//When eomploying OpenGL run using: g++ -std=c++17 main.cpp gpu_renderer.cpp simulation.cpp -o main -lGLEW -lGL -lGLU -lsfml-graphics -lsfml-window -lsfml-system
#include "gpu_renderer.h"
#include "simulation.h"

int main() {
    constexpr int particleCount = 20000;
    auto particles = createParticles(particleCount);
    GPURenderer renderer;
    return renderer.run(particles);
}
