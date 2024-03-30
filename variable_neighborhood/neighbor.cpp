#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "leasing.cpp"

struct Neighbor {
    long long fitness=-1;
    std::vector<std::set<std::pair<int, int>>> dna;
    std::vector<std::set<int>> current_facilities;
    std::vector<std::pair<int, int>> index_best;
    std::vector<long long> clients_best;
    LeasingProblem problem;

    Neighbor(LeasingProblem const &_problem) {
        problem = _problem;
        dna = std::vector<std::set<std::pair<int, int>>>(problem.T);
        current_facilities = std::vector<std::set<int>>(problem.T);
        index_best = std::vector<std::pair<int, int>>(problem.V);
        clients_best = std::vector<long long>(problem.V);
    }

    Neighbor(LeasingProblem const &_problem, std::vector<std::set<std::pair<int, int>>> _dna, long long _fitness) {
        problem = _problem;
        index_best = std::vector<std::pair<int, int>>(problem.V);
        clients_best = std::vector<long long>(problem.V);

        dna = _dna;
        current_facilities = std::vector<std::set<int>>(problem.T);
        fitness = _fitness;
        get_facilities();
    }

    void shake() {

    }

    void local_search() {
        std::vector<std::set<std::pair<int, int>>> new_sol(problem.T);

        int cnt = 0;
        for (int t = 0; t < problem.T; t++) {
            cnt += (int)dna[t].size();
        }
        for (int i = 0; i < cnt; i++) index_best[i] = {0, i};

        for (int t = 0; t < problem.T; t++) {
            for (int v = 0; v < problem.V; v++) clients_best[v] = INF;
            for (int client: problem.clients[t]) {
                int best = 0;
                for (int i = 0, c = 0; i < problem.T; i++) {
                    if (dna[i].empty()) {
                        c++; continue;
                    }
                    for (auto [u, l]: dna[i]) {
                        int nd = std::min(problem.T, i + problem.center_types[l]);
                        if (i <= t && t < nd) {
                            if (problem.graph[client][u] < clients_best[client]) {
                                clients_best[client] = problem.graph[client][u];
                                best = c;
                            }
                        }
                        c++;
                    }
                }
                index_best[best].first++;
            }
        }

        std::sort(index_best.begin(), index_best.end());
        int choosen = rand_i() % cnt; 

        for (int i = 0; i < choosen; i++) {
            int id = index_best[i].second;
            for (int t = 0, c = 0; t < problem.T; t++) {
                if (dna[t].empty()) continue;
                if (c + (int)dna.size() <= id) {
                    for (auto [u, l]: dna[t]) {
                        if (c == id) {
                            new_sol[t].insert({u, l});
                        }
                        c++;
                    }  
                    break;             
                }
                c += (int)dna.size();
            }
        }

        for (int t = 0; t < problem.T; t++) {
            current_facilities[t].clear();
        }
        bool ok = true;
        int t = 0;
        while (ok) {
            int l = rand_i() % problem.L;
            ok = false;
            for (int i = t; i < problem.T; i++) {
                int st = i;
                int nd = std::min(problem.T, st + problem.center_types[l]);
                int mx = 0;
                for (int j = st; j < nd; j++) {
                    mx = std::max(mx, (int)current_facilities[st].size());
                }
                if (mx < problem.K) {
                    std::set<int> available;
                    for (int v = 0; v < problem.V; v++) {
                        bool has = false;
                        for (int j = st; j < nd && !has; j++){
                            if (current_facilities[j].find(v) != current_facilities[j].end()) {
                                has = true;
                            }
                        }
                        if (!has) {
                            available.insert(v);
                        }
                    }
                    if (available.empty()){
                        continue;
                    }
                    int op = rand_i() % (int)available.size(), v;
                    for (int j: available) {
                        if (!op) {
                            v = j;
                        }
                        op--;
                    }
                    for (int u = st; u < nd; u++){
                        current_facilities[u].insert(v);
                    }
                    new_sol[i].insert({v, l});
                    ok = true;
                    t = i+1;
                    break;
                }
                i = nd;
            }
        }

        for (t = 0; t < problem.T; t++) {
            dna[t] = new_sol[t];
        }

        get_fitness();
    }

    void get_facilities() {
        for (int t = 0; t < problem.T; t++) {
            for (auto [v, l]: dna[t]) {
                for (int i = t; i < std::min(problem.T, t + problem.center_types[l]); i++) {
                    current_facilities[i].insert(v);
                }
            }
        }
    }

    int get_fitness() {
        fitness = 0;

        get_facilities();

        for (int t = 0; t < problem.T; t++) {
            for (int client: problem.clients[t]) {
                long long mn = INF;
                for (int facility: current_facilities[t]) {
                    if (facility == -1) break;
                    mn = std::min(mn, problem.graph[client][facility]);
                }
                fitness += mn;
            }
        }
        
        return fitness;
    }

    friend bool operator < (Neighbor const &lhs, Neighbor const &rhs){
        return lhs.fitness < rhs.fitness;
    }

    friend std::ostream &operator <<(std::ostream &os, const Neighbor &sol) {
        os << "Fitness: " << sol.fitness;
        return os;
    }
};