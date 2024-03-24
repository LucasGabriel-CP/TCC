#pragma once

#include <vector>
#include <random>
#include <algorithm>

#include "Cable.cpp"
#include "Point.cpp"

struct Turbine {
    int id, total_prod;
    Point pos;
    Cable cable;
    std::vector<int> connetions;
    Turbine() : id(-1), total_prod(0), connetions({}), cable() { }
    Turbine(int _id, int _total_prod, Point _pos, std::vector<int> cons = {})
        : id(_id), total_prod(_total_prod), pos(_pos), connetions(cons) { }

    Turbine &operator+=(Turbine rhs) & {
        total_prod += rhs.total_prod;
        connetions.push_back(rhs.id);
        return *this;
    }

    friend Turbine operator+(Turbine const &lhs, Turbine const &rhs) {
        Turbine res = lhs;
        res += rhs;
        return res;
    }

    bool check(Turbine const &rhs) {
        return (total_prod + rhs.total_prod) <= cable.capacity;
    }
};
