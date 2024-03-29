#pragma once

#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <map>
#include "util/DebugTemplate.cpp"

int V, T, L, K, E;
std::vector<std::vector<std::pair<int, int>>> graph;      // dist(i, j)
std::vector<std::vector<int>> clients;          // Clientes no instante t
std::vector<int> center_types;                  // Tipos de centro

void read(std::string const &instance_path) {
    std::ifstream in;

    debug(instance_path);
    in.open(instance_path);

    in >> V >> T >> L >> K;

    graph.assign(V, std::vector<std::pair<int, int>>(V));
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            graph[i][j].first = j;
            in >> graph[i][j].second;
        }
        std::sort(graph[i].begin(), graph[i].end(), [](auto a, auto b){
            return a.second < b.second;
        });
    }

    center_types.resize(L);
    for (int i = 0; i < L; i++) {
        in >> center_types[i];
    }
    debug(center_types);

    clients.resize(T);
    for (int t = 0; t < T; t++) {
        int dt; in >> dt;
        clients[t].resize(dt);
        for (int j = 0; j < dt; j++) {
            in >> clients[t][j];
        }
        debug(clients[t]);
    }

    in.close();
}

