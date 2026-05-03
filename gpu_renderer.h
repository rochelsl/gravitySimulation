#pragma once

#include <vector>
#include "particle.h"

class GPURenderer {
public:
    int run(const std::vector<Particle>& particles);
};
