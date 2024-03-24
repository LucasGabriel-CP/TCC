#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include "Cable.cpp"
#include "Point.cpp"
#include "DebugTemplate.cpp"
#include "Constants.cpp"

struct Individuo {
    Point substation;
    int n_nodes, n_steiner;
    std::vector<std::vector<int>> branches;
    std::vector<int> in_deg, out_deg;
    float fitness;

    Individuo() { }

    Individuo(Point _substation, int _n_nodes) {
        substation = _substation;
        n_nodes = _n_nodes;
        n_steiner = 0;
        branches = std::vector<std::vector<int>>(n_nodes);
        in_deg = std::vector<int>(n_nodes);
        out_deg = std::vector<int>(n_nodes);
    }

    Individuo(
        Point _substation, int _n_nodes, std::vector<std::vector<int>> _branches,
        std::vector<int> _in_deg, std::vector<int> _out_deg, int _n_steiner
    ) {
        substation = _substation;
        n_nodes = _n_nodes;
        n_steiner = _n_steiner;
        branches = _branches;
        in_deg = _in_deg;
        out_deg = _out_deg;
    }

    void build(int n_branches = 0) {
        if (!n_branches) n_branches = n_nodes;
        std::vector<int> roots(n_nodes);
        std::iota(roots.begin(), roots.end(), 0);

        std::shuffle(roots.begin(), roots.end(), g);

        roots = std::vector<int>(roots.begin(), roots.begin() + n_branches - 1);
        debug(roots);

        branches[0] = roots;
        
        int avaliable = n_nodes;
    }

    void mutate() {

        fitness = get_fitness();
    }

    static std::pair<Individuo, Individuo> cross(Individuo const &p1, Individuo const &p2) {
        return {Individuo(), Individuo()};
    }

    float get_fitness() {
        return 0.0;
    }
    friend bool operator < (const Individuo& lhs, const Individuo& rhs){
        return lhs.fitness < rhs.fitness;
    }
};
