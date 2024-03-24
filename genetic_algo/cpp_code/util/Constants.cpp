#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <mutex>

const double EPS = 1e-6;
const float FLOAT_INF = std::numeric_limits<float>::max();
std::mutex mtx;

#define rand() uid(g)
std::mt19937 g(std::chrono::high_resolution_clock::now().time_since_epoch().count());
std::uniform_real_distribution<float> uid(0.0, 1.0);

#define NIL -1