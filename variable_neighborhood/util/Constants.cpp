#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <mutex>

std::mutex mtx;

#define rand() uid(g)
thread_local std::mt19937 g(std::chrono::high_resolution_clock::now().time_since_epoch().count());
thread_local std::uniform_real_distribution<double> uid(0.0, 1.0);

#define rand_i() uid_i(g)
thread_local std::uniform_int_distribution<int> uid_i(0, std::numeric_limits<int>::max());             // range

const double EPS = 1e-6;
const double DOUBLE_INF = std::numeric_limits<double>::max();

const long long INF = std::numeric_limits<long long>::max();
const long long PENALTY = 2;

int f_cmp(double const &lhs, double const &rhs) {
    if (fabs(lhs - rhs) < EPS) return 0;
    return lhs < rhs ? -1 : 1;
}