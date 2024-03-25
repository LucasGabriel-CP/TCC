#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include "util/Cable.cpp"
#include "util/Point.cpp"
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "util/Turbine.cpp"

std::vector<std::vector<double>> graph;
std::vector<double> branches;
std::vector<Cable> cables;
std::vector<Turbine> turbinies;
Point substation;
int N, T, C;

double dist(int i, int j) {
    if (i == NIL) return branches[j];
    if (j == NIL) return branches[i];
    return graph[i][j];
}

void init_graph() {
    graph.assign(N, std::vector<double>(N, 0.0));
    branches.assign(N, 0.0);

    for (int i = 0; i < N; i++) {
        branches[i] = dist(substation, turbinies[i].pos);
        for (int j = i + 1; j < N; j++) {
            graph[i][j] = graph[j][i] = dist(turbinies[i].pos, turbinies[j].pos);
        }
    }
}

Cable best_cable(int energy) {
    for (Cable &cable: cables) {
        if (cable.capacity >= energy) {
            return cable;
        }
    }
    return Cable(T, N, ENERGY_OVERFLOW_PENALTY, 999);
    // assert(false && "Energy overflow");
}

int best_cable_price(int energy) {
    return best_cable(energy).price;
}

void read(std::string const &inf_file, std::string const &turb_file, std::string const &cable_file) {
    std::ifstream in;

    debug(inf_file);
    in.open(inf_file);
    in >> N >> C;
    in.close();

    debug(turb_file);
    in.open(turb_file);

    in >> substation;
    turbinies.resize(N);
    int garb; in >> garb;
    for (int i = 0; i < N; i++) {
        in >> turbinies[i].pos >> turbinies[i].total_prod;
        turbinies[i].id = i;
        turbinies[i].total_prod = 1;
    }

    in.close();

    debug(cable_file);
    in.open(cable_file);
    assert(in.good());

    in >> T;
    cables.resize(T);
    for (int i = 0; i < T; i++) {
        in >> cables[i].capacity >> cables[i].price >> cables[i].max_usage;
        cables[i].id = i;
    }

    in.close();
}