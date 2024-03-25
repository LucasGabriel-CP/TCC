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

#define rand_i() uid_i(g)
std::uniform_int_distribution<int> uid_i(0, std::numeric_limits<int>::max());             // range

int f_cmp(float const &lhs, float const &rhs) {
    if (fabs(lhs - rhs) < EPS) return 0;
    return lhs < rhs ? -1 : 1;
}