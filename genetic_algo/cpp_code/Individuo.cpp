#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include "util/Cable.cpp"
#include "util/Point.cpp"
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "util/Turbine.cpp"
#include "SWFCR.cpp"

struct Individuo {
    Point substation;
    int n_nodes, n_steiner, sub_connections;
    std::vector<std::vector<int>> branches;
    std::vector<int> in_deg, out_deg, go;
    std::vector<Cable> cables;
    float fitness;

    Individuo() { }

    Individuo(Point _substation, int _n_nodes) {
        substation = _substation;
        n_nodes = _n_nodes;
        n_steiner = 0;
        branches = std::vector<std::vector<int>>(n_nodes);
        in_deg = std::vector<int>(n_nodes);
        out_deg = std::vector<int>(n_nodes);
        go = std::vector<int>(n_nodes);
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
        go = std::vector<int>(n_nodes);
    }

    std::vector<std::vector<int>> create_layers(std::vector<Turbine> const &turbinies, std::vector<int> const& order, int start) {
        std::vector<int> layer;
        std::vector<std::vector<int>> layers;
        int total_energy = 0;
        int biggest_capacity = cables.back().capacity;
        for (int i = start - n_nodes; ;i++) {
            int at = order[i < 0 ? i + n_nodes : i];
            if (i == start || total_energy + turbinies[at].total_prod > biggest_capacity) {
                layers.push_back(layer);
                total_energy = 0;
            }
            if (i == start) return layers;
            total_energy += turbinies[at].total_prod;
        }
        assert(false);
    }

    void build(std::vector<Turbine> const &turbinies, std::vector<int> const& order, int start) {
        std::vector<std::vector<int>> layers = create_layers(turbinies, order, start);
        for (std::vector<int> layer: layers) {
            debug(layer);
            int u = layer.back();
            for (int i = (int)layer.size() - 2; i >= 0; i--) {
                go[u] = layer[i];
                u = layer[i];
            }
            go[u] = NIL;
            sub_connections++;
        }
        get_fitness();
    }

    void mutate() {
        
        get_fitness();
    }

    static std::pair<Individuo, Individuo> cross(Individuo const &p1, Individuo const &p2) {

        return {Individuo(), Individuo()};
    }

    float get_fitness() {
        fitness = 0.0;

        for (int i = 0; i < n_nodes; i++) {
            fitness += dist(i, go[i]) * best_cable_price(turbinies[i].total_prod);
        }
        
        return fitness;
    }
    friend bool operator < (const Individuo& lhs, const Individuo& rhs){
        return lhs.fitness < rhs.fitness;
    }
};
