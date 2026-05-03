# Gravity Simulation

A simple 2D particle-based gravity simulation written in C++ using SFML, OpenGL, GLEW, and GPU compute shaders. The program initializes thousands of particles, computes their mutual gravitational motion approximately on a density grid, and renders them as colored points in real time.

## Capabilities

- Simulates many particles; default: `20000`
- Uses GPU compute shaders for density deposition, force-field calculation, force sampling, and integration
- Supports periodic boundary conditions: particles leaving one side of the window re-enter from the opposite side
- Initializes particles either from a normal distribution or from a jittered grid
- Adds a small initial rotational velocity field
- Colors particles by speed: slow particles appear blue, faster particles become white
- Uses a softened gravitational potential to avoid singular forces at very small distances

## Physics / Mathematics

Each particle has position, velocity, acceleration, mass, and radius:

$\mathbf{x}_i,\quad \mathbf{v}_i,\quad \mathbf{a}_i,\quad m_i,\quad r_i$

The simulation approximates gravity by first assigning particles to a 2D density grid. Each grid cell acts like a lumped mass:

$M_j = N_j m$

where $N_j$ is the number of particles in cell $j$, and $m$ is the particle mass.

The gravitational acceleration from a cell is computed as:

$$
\mathbf{a}_i =
G M_j
\frac{\mathbf{r}_{ij}}
{\left(|\mathbf{r}_{ij}|^2 + \epsilon^2\right)^{3/2}}
$$

where:

- $G$ is the gravitational constant
- $\mathbf{r}_{ij}$ is the displacement from particle or cell $i$ to cell $j$
- $\epsilon$ is a softening length that prevents infinite acceleration at zero distance

The potential energy used during initialization is approximately:

$U =-\sum_{i<j} \frac{G m_i m_j}{\sqrt{|\mathbf{r}_{ij}|^2 + 25}}$

The kinetic energy is:

$$
T =
\sum_i
\frac{1}{2} m_i |\mathbf{v}_i|^2
$$

Particle motion is advanced with a simple explicit update:

$\mathbf{v}_{t+\Delta t}=\mathbf{v}_t + \mathbf{a}_t \Delta t$

$$\mathbf{x}_{t+\Delta t}=\mathbf{x}_t + \mathbf{v}_{t+\Delta t} \Delta t$$

A small damping factor is applied each step:

$\mathbf{v} \leftarrow 0.99999 \mathbf{v}$

## Important Parameters

Most simulation parameters are defined in `constants.h`.

| Parameter | Meaning | Effect |
|---|---|---|
| `kWidth`, `kHeight` | Simulation window size | Larger values increase the spatial domain |
| `particleCount` | Number of simulated particles | More particles create richer structures but cost more GPU work |
| `kGConst` | Gravitational constant | Higher values make attraction stronger and collapse faster |
| `kParticleMass` | Mass of each particle | Higher mass increases gravitational attraction |
| `kDt` | Time step | Larger values speed up motion but can reduce stability |
| `kInitType` | Initial particle layout | `"normal"` gives a clustered random cloud; `"grid"` gives a jittered lattice |
| `kGridW`, `kGridH` | Density-grid resolution | Higher values give a finer force approximation but require more computation |
| `kComputeLocalSize` | Compute shader workgroup size | Affects GPU execution layout |
| `eps` | Gravity softening length | Larger values smooth the force field and reduce close-range collapse |
| `damping` | Velocity damping | Lower values remove energy faster |
| `rotationStrength` | Initial swirl strength | Higher values add more orbital/rotational motion |
| `maxSpeed` | Rendering color scale | Lower values make particles appear bright at lower speeds |
| `uPointSize` | Rendered particle size | Larger values make particles easier to see |

## Build
Compile with:

```bash
g++ -std=c++17 main.cpp gpu_renderer.cpp simulation.cpp -o main \
  -lGLEW -lGL -lGLU -lsfml-graphics -lsfml-window -lsfml-system
