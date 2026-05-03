#include "gpu_renderer.h"

#include <GL/glew.h>
#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>

#include "constants.h"

namespace {
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
    glDeleteShader(vs); glDeleteShader(fs);
    return program;
}

GLuint createComputeProgram(const char* computeSrc) {
    GLuint cs = compileShader(GL_COMPUTE_SHADER, computeSrc);
    GLuint program = glCreateProgram();
    glAttachShader(program, cs);
    glLinkProgram(program);
    glDeleteShader(cs);
    return program;
}
}

int GPURenderer::run(const std::vector<Particle>& particles) {
    sf::ContextSettings settings; settings.majorVersion = 4; settings.minorVersion = 3;
    auto window = sf::RenderWindow(sf::VideoMode(kWidth, kHeight), "Particle Simulation", sf::Style::Default, settings);
    window.setActive(true);
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) return -1;

    const char* vertexShaderSrc = R"(#version 430 core
layout(std430, binding = 0) buffer Positions { vec4 posMass[]; };
layout(std430, binding = 1) buffer Velocities { vec4 velRadius[]; };
uniform float uWidth; uniform float uHeight; uniform float uPointSize; out float vSpeed;
void main(){ vec2 pos = posMass[gl_VertexID].xy; vec2 vel = velRadius[gl_VertexID].xy; vSpeed=length(vel);
float x=(pos.x/uWidth)*2.0-1.0; float y=1.0-(pos.y/uHeight)*2.0; gl_Position=vec4(x,y,0,1); gl_PointSize=uPointSize; })";
    const char* fragmentShaderSrc = R"(#version 430 core
in float vSpeed; out vec4 FragColor; uniform float maxSpeed;
void main(){ float t=clamp(vSpeed/maxSpeed,0.0,1.0); vec3 color=mix(vec3(0.2,0.4,1.0),vec3(1.0),t); FragColor=vec4(color,1.0);} )";
    const char* integrateShaderSrc = R"(#version 430 core
layout(local_size_x = 256) in; layout(std430, binding = 0) buffer Positions { vec4 posMass[]; }; layout(std430, binding = 1) buffer Velocities { vec4 velRadius[]; }; layout(std430, binding = 2) buffer Accelerations { vec4 acc[]; };
uniform float dt; uniform float width; uniform float height; uniform float damping; uniform uint numParticles;
void main(){ uint i=gl_GlobalInvocationID.x; if(i>=numParticles) return; vec2 pos=posMass[i].xy; vec2 vel=velRadius[i].xy; vec2 a=acc[i].xy; vel+=a*dt; vel*=damping; pos+=vel*dt; if(pos.x<0.0) pos.x+=width; if(pos.x>=width) pos.x-=width; if(pos.y<0.0) pos.y+=height; if(pos.y>=height) pos.y-=height; posMass[i].xy=pos; velRadius[i].xy=vel; })";
    const char* clearGridShaderSrc = R"(#version 430 core
layout(local_size_x = 256) in; layout(std430, binding = 3) buffer DensityGrid { uint density[]; }; uniform uint gridSize;
void main(){ uint i=gl_GlobalInvocationID.x; if(i>=gridSize) return; density[i]=0u;})";
    const char* depositShaderSrc = R"(#version 430 core
layout(local_size_x = 256) in; layout(std430, binding = 0) buffer Positions { vec4 posMass[]; }; layout(std430, binding = 3) buffer DensityGrid { uint density[]; };
uniform float width; uniform float height; uniform uint gridW; uniform uint gridH; uniform uint numParticles;
void main(){ uint i=gl_GlobalInvocationID.x; if(i>=numParticles) return; vec2 pos=posMass[i].xy; uint cellX=uint(clamp(pos.x/width*float(gridW),0.0,float(gridW-1))); uint cellY=uint(clamp(pos.y/height*float(gridH),0.0,float(gridH-1))); atomicAdd(density[cellY*gridW+cellX],1u);} )";
    const char* computeGridForceShaderSrc = R"(#version 430 core
layout(local_size_x = 256) in; layout(std430, binding = 3) buffer DensityGrid { uint density[]; }; layout(std430, binding = 4) buffer ForceGrid { vec4 force[]; };
uniform uint gridW; uniform uint gridH; uniform float width; uniform float height; uniform float gConst; uniform float particleMass; uniform float eps;
void main(){ uint idx=gl_GlobalInvocationID.x; uint gridSize=gridW*gridH; if(idx>=gridSize) return; uint x=idx%gridW; uint y=idx/gridW; vec2 cellSize=vec2(width/float(gridW),height/float(gridH)); vec2 posA=vec2((float(x)+0.5)*cellSize.x,(float(y)+0.5)*cellSize.y); vec2 acc=vec2(0.0); for(uint j=0;j<gridSize;++j){ uint count=density[j]; if(count==0u||j==idx) continue; uint xj=j%gridW; uint yj=j/gridW; vec2 posB=vec2((float(xj)+0.5)*cellSize.x,(float(yj)+0.5)*cellSize.y); vec2 r=posB-posA; if(r.x>width*0.5) r.x-=width; if(r.x<-width*0.5) r.x+=width; if(r.y>height*0.5) r.y-=height; if(r.y<-height*0.5) r.y+=height; float mass=float(count)*particleMass; float dist2=dot(r,r)+eps*eps; float invDist=inversesqrt(dist2); float invDist3=invDist*invDist*invDist; acc+=gConst*mass*r*invDist3;} force[idx]=vec4(acc,0,0);} )";
    const char* sampleForceShaderSrc = R"(#version 430 core
layout(local_size_x = 256) in; layout(std430, binding = 0) buffer Positions { vec4 posMass[]; }; layout(std430, binding = 2) buffer Accelerations { vec4 acc[]; }; layout(std430, binding = 4) buffer ForceGrid { vec4 force[]; };
uniform float width; uniform float height; uniform uint gridW; uniform uint gridH; uniform uint numParticles;
void main(){ uint i=gl_GlobalInvocationID.x; if(i>=numParticles) return; vec2 pos=posMass[i].xy; uint cellX=uint(clamp(pos.x/width*float(gridW),0.0,float(gridW-1))); uint cellY=uint(clamp(pos.y/height*float(gridH),0.0,float(gridH-1))); acc[i]=force[cellY*gridW+cellX]; })";

    GLuint renderProgram = createProgram(vertexShaderSrc, fragmentShaderSrc);
    GLuint integrateProgram = createComputeProgram(integrateShaderSrc);
    GLuint clearGridProgram = createComputeProgram(clearGridShaderSrc);
    GLuint depositProgram = createComputeProgram(depositShaderSrc);
    GLuint computeGridForceProgram = createComputeProgram(computeGridForceShaderSrc);
    GLuint sampleForceProgram = createComputeProgram(sampleForceShaderSrc);

    std::vector<GPUPosition> gpuPositions(particles.size());
    std::vector<GPUVelocity> gpuVelocities(particles.size());
    std::vector<GPUAcceleration> gpuAccelerations(particles.size());
    for (size_t i=0;i<particles.size();++i){ gpuPositions[i]={particles[i].position.x,particles[i].position.y,0,particles[i].mass}; gpuVelocities[i]={particles[i].velocity.x,particles[i].velocity.y,0,particles[i].radius}; gpuAccelerations[i]={0,0,0,0}; }

    GLuint posSSBO, velSSBO, accSSBO, densitySSBO, forceGridSSBO, vao;
    glGenVertexArrays(1, &vao); glBindVertexArray(vao);
    glGenBuffers(1, &posSSBO); glBindBuffer(GL_SHADER_STORAGE_BUFFER, posSSBO); glBufferData(GL_SHADER_STORAGE_BUFFER, gpuPositions.size()*sizeof(GPUPosition), gpuPositions.data(), GL_DYNAMIC_DRAW); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,posSSBO);
    glGenBuffers(1, &velSSBO); glBindBuffer(GL_SHADER_STORAGE_BUFFER, velSSBO); glBufferData(GL_SHADER_STORAGE_BUFFER, gpuVelocities.size()*sizeof(GPUVelocity), gpuVelocities.data(), GL_DYNAMIC_DRAW); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,velSSBO);
    glGenBuffers(1, &accSSBO); glBindBuffer(GL_SHADER_STORAGE_BUFFER, accSSBO); glBufferData(GL_SHADER_STORAGE_BUFFER, gpuAccelerations.size()*sizeof(GPUAcceleration), gpuAccelerations.data(), GL_DYNAMIC_DRAW); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,2,accSSBO);
    std::vector<unsigned int> densityGrid(kGridSize,0); glGenBuffers(1,&densitySSBO); glBindBuffer(GL_SHADER_STORAGE_BUFFER,densitySSBO); glBufferData(GL_SHADER_STORAGE_BUFFER,densityGrid.size()*sizeof(unsigned int),densityGrid.data(),GL_DYNAMIC_DRAW); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,3,densitySSBO);
    std::vector<GPUForce> forceGrid(kGridSize); glGenBuffers(1,&forceGridSSBO); glBindBuffer(GL_SHADER_STORAGE_BUFFER,forceGridSSBO); glBufferData(GL_SHADER_STORAGE_BUFFER,forceGrid.size()*sizeof(GPUForce),forceGrid.data(),GL_DYNAMIC_DRAW); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,4,forceGridSSBO);

    while(window.isOpen()){
        sf::Event event; while(window.pollEvent(event)){ if(event.type==sf::Event::Closed) window.close(); }
        GLuint particleGroups = static_cast<GLuint>((particles.size()+255)/256); GLuint gridGroups = static_cast<GLuint>((kGridSize+255)/256);
        glUseProgram(clearGridProgram); glUniform1ui(glGetUniformLocation(clearGridProgram,"gridSize"), static_cast<unsigned int>(kGridSize)); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,3,densitySSBO); glDispatchCompute(gridGroups,1,1); glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUseProgram(depositProgram); glUniform1f(glGetUniformLocation(depositProgram,"width"), static_cast<float>(kWidth)); glUniform1f(glGetUniformLocation(depositProgram,"height"), static_cast<float>(kHeight)); glUniform1ui(glGetUniformLocation(depositProgram,"gridW"),kGridW); glUniform1ui(glGetUniformLocation(depositProgram,"gridH"),kGridH); glUniform1ui(glGetUniformLocation(depositProgram,"numParticles"), static_cast<unsigned int>(particles.size())); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,posSSBO); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,3,densitySSBO); glDispatchCompute(particleGroups,1,1); glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUseProgram(computeGridForceProgram); glUniform1ui(glGetUniformLocation(computeGridForceProgram,"gridW"),kGridW); glUniform1ui(glGetUniformLocation(computeGridForceProgram,"gridH"),kGridH); glUniform1f(glGetUniformLocation(computeGridForceProgram,"width"), static_cast<float>(kWidth)); glUniform1f(glGetUniformLocation(computeGridForceProgram,"height"), static_cast<float>(kHeight)); glUniform1f(glGetUniformLocation(computeGridForceProgram,"gConst"),kGConst); glUniform1f(glGetUniformLocation(computeGridForceProgram,"particleMass"),kParticleMass); glUniform1f(glGetUniformLocation(computeGridForceProgram,"eps"),20.0f); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,3,densitySSBO); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,4,forceGridSSBO); glDispatchCompute(gridGroups,1,1); glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUseProgram(sampleForceProgram); glUniform1f(glGetUniformLocation(sampleForceProgram,"width"),static_cast<float>(kWidth)); glUniform1f(glGetUniformLocation(sampleForceProgram,"height"),static_cast<float>(kHeight)); glUniform1ui(glGetUniformLocation(sampleForceProgram,"gridW"),kGridW); glUniform1ui(glGetUniformLocation(sampleForceProgram,"gridH"),kGridH); glUniform1ui(glGetUniformLocation(sampleForceProgram,"numParticles"),static_cast<unsigned int>(particles.size())); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,posSSBO); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,2,accSSBO); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,4,forceGridSSBO); glDispatchCompute(particleGroups,1,1); glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUseProgram(integrateProgram); glUniform1f(glGetUniformLocation(integrateProgram,"dt"),kDt); glUniform1f(glGetUniformLocation(integrateProgram,"width"),static_cast<float>(kWidth)); glUniform1f(glGetUniformLocation(integrateProgram,"height"),static_cast<float>(kHeight)); glUniform1f(glGetUniformLocation(integrateProgram,"damping"),0.99999f); glUniform1ui(glGetUniformLocation(integrateProgram,"numParticles"), static_cast<unsigned int>(particles.size())); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,posSSBO); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,velSSBO); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,2,accSSBO); glDispatchCompute(particleGroups,1,1); glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT);
        glViewport(0,0,kWidth,kHeight); glClearColor(0,0,0,1); glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(renderProgram); glUniform1f(glGetUniformLocation(renderProgram,"uWidth"),static_cast<float>(kWidth)); glUniform1f(glGetUniformLocation(renderProgram,"uHeight"),static_cast<float>(kHeight)); glUniform1f(glGetUniformLocation(renderProgram,"uPointSize"),1.5f); glUniform1f(glGetUniformLocation(renderProgram,"maxSpeed"),500.0f); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,posSSBO); glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,velSSBO); glBindVertexArray(vao); glDrawArrays(GL_POINTS,0,static_cast<GLsizei>(particles.size()));
        window.display();
    }
    return 0;
}
