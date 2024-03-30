#pragma once

#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <map>
#include "util/DebugTemplate.cpp"


struct LeasingProblem {
    int V, T, L, K;
    std::vector<std::vector<long long>> graph;
    std::vector<std::vector<std::pair<int, long long>>> adj;
    std::vector<std::vector<int>> clients;          // Clientes no instante t
    std::vector<int> center_types;                  // Tipos de centro
    std::string version;

    double mean_l;

    LeasingProblem() { }

    LeasingProblem(std::string const &instance_path, std::string const &_version) {
        read(instance_path);
        version = _version;
    }

    void read(std::string const &instance_path) {
        std::ifstream in;

        in.open(instance_path);

        in >> V >> T >> L >> K;

        graph.assign(V, std::vector<long long>(V));
        adj.resize(V);
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                in >> graph[i][j];
                adj[i].push_back({j, graph[i][j]});
            }
            std::sort(adj[i].begin(), adj[i].end(), [](auto a, auto b){
                return a.second < b.second;
            });
        }

        center_types.resize(L);
        mean_l = 0;
        for (int i = 0; i < L; i++) {
            in >> center_types[i];
            mean_l += i;
        }
        mean_l /= L;

        clients.resize(T);
        for (int t = 0; t < T; t++) {
            int dt; in >> dt;
            clients[t].resize(dt);
            for (int j = 0; j < dt; j++) {
                in >> clients[t][j];
            }
        }

        in.close();
    }    
};