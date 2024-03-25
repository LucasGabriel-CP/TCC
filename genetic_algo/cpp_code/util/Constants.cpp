#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <mutex>

#define rand() uid(g)
std::mt19937 g(std::chrono::high_resolution_clock::now().time_since_epoch().count());
std::uniform_real_distribution<double> uid(0.0, 1.0);

#define rand_i() uid_i(g)
std::uniform_int_distribution<int> uid_i(0, std::numeric_limits<int>::max());             // range

#define NIL -1
#define CLICE_PENALTY 1e3
#define NOT_REACH_PENALTY 1e10
#define C_OVERFLOW_PENALTY 1e6
#define DISCONNECTED_PENALTY 1e8
#define CROSS_PENLATY 1e3
#define ENERGY_OVERFLOW_PENALTY 1e8

const double EPS = 1e-6;
const double DOUBLE_INF = std::numeric_limits<double>::max();

int f_cmp(double const &lhs, double const &rhs) {
    if (fabs(lhs - rhs) < EPS) return 0;
    return lhs < rhs ? -1 : 1;
}