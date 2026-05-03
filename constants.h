#pragma once

#include <string>

inline constexpr int kWidth = 1500;
inline constexpr int kHeight = 1500;

inline constexpr float kGConst = 10.0f;
inline constexpr float kParticleMass = 5.0f;

inline constexpr char kInitType[] = "normal";

inline constexpr int kGridW = 64;
inline constexpr int kGridH = 64;
inline constexpr int kGridSize = kGridW * kGridH;

inline constexpr int kComputeLocalSize = 256;
inline constexpr float kDt = 0.0005f;
