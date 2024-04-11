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

#define MULTI_FACILITY_PENALTY 500000LL
#define NO_FACILITY_PENALTY 5000LL
#define PK_FACILITIES_PENALTY 5000LL
#define ALL_PENALTIES (MULTI_FACILITY_PENALTY + PK_FACILITIES_PENALTY * PK_FACILITIES_PENALTY + NO_FACILITY_PENALTY * NO_FACILITY_PENALTY * NO_FACILITY_PENALTY)


int f_cmp(double const &lhs, double const &rhs) {
    if (fabs(lhs - rhs) < EPS) return 0;
    return lhs < rhs ? -1 : 1;
}

double get_rand(double l, double r) {
    std::uniform_real_distribution<double> distri(l, r);
    return distri(g);
}