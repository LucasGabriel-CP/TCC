#pragma once

#include <vector>
#include <set>
#include "util/DebugTemplate.cpp"
#include "util/segtree.cpp"
#include "util/Constants.cpp"
#include "leasing.cpp"


struct InitialSolution {
    long long fitness=-1;
    std::vector<std::set<std::pair<int, int>>> dna;
    std::vector<std::set<int>> current_facilities;
    std::vector<std::pair<int, int>> rec;
    std::vector<std::vector<bool>> vis;
    std::vector<int> aproximation_server, dp;
    std::vector<SegTree> facility_seg;
    LeasingProblem problem;


    InitialSolution(LeasingProblem const &_problem) {
        problem = _problem;

        dna = std::vector<std::set<std::pair<int, int>>>(problem.T);
        current_facilities = std::vector<std::set<int>>(problem.T);
        rec = std::vector<std::pair<int, int>>(problem.T);
        aproximation_server = std::vector<int>(problem.K, -1);
        dp = std::vector<int>(problem.T, -1);
        vis = std::vector<std::vector<bool>>(problem.T, std::vector<bool>(problem.V, false));
        facility_seg = std::vector<SegTree>(problem.V, SegTree(problem.T));
    }

    void dorit(int t) {
        int low = 0, high = problem.V-1;
        std::vector<int> ans;
        while (low < high) {
            int mid = (low + high) / 2;
            std::vector<int> S;
            std::set<int> set_V;
            for (int i = 0; i < problem.V; i++) set_V.insert(i);
            while (!set_V.empty()) {
                int x = *set_V.begin();
                S.push_back(x);
                for (int i = 0; i <= mid; i++) {
                    int v = problem.adj[x][i].first;
                    if (set_V.find(v) != set_V.end()) {
                        set_V.erase(v);
                        for (int j = 0; j <= mid; j++) {
                            int z = problem.adj[v][j].first;
                            if (set_V.find(z) != set_V.end()) {
                                set_V.erase(z);
                            }
                        }
                    }
                }
                if ((int)S.size() <= problem.K) {
                    ans = S;
                    high = mid;
                }
                else {
                    low = mid + 1;
                }
            }
        }
        for (int i: ans) {
            facility_seg[i].update(t, t, 1);
        }
    }

    void gon(int t) {
        aproximation_server[0] = rand_i() % problem.V;
        int k = 1;
        while(k < problem.K) {
            int next_server = -1, distance = 0;
            for(int v = 0; v < problem.V; v++) {
                long long mn = INF;
                bool is_in_the_set = false;
                for(int x: aproximation_server){
                    if (x == -1) break;
                    if (v == x) {
                        is_in_the_set = true;
                        break;
                    }
                    mn = std::min(mn, problem.graph[v][x]);
                }
                if(mn > distance && !is_in_the_set) {
                    next_server = v;
                    distance = mn;
                }
            }
            aproximation_server[k++] = next_server;
        }
        for (int i: aproximation_server) {
            facility_seg[i].update(t, t, 1);
        }
    }

    int jorge_algos(int t) {
        if (t >= problem.T) return 0;
        auto &d = dp[t];
        if (d != -1) return d;
        d = 0;
        for (int l = 0; l < problem.L; l++) {
            int lim = std::min(problem.T, t + problem.center_types[l]);
            for (int v = 0; v < problem.V; v++) {
                if (vis[t][v]) continue;
                int qnt = facility_seg[v].query(t, lim-1);
                int aux = jorge_algos(lim) + qnt;
                if (aux > d)  {
                    rec[t] = {l, v};
                    d = aux;
                }
            }
        }
        return d;
    }


    void build(bool use_dorit = false) {
        for (int t = 0; t < problem.T; t++) {
            aproximation_server = std::vector<int>(problem.K, -1);
            if (use_dorit) {
                dorit(t);
            }
            else {
                gon(t);
            }
        }
        for (int k = 0; k < problem.K; k++) {
            for (int &t: dp) t = -1;
            jorge_algos(0);
            int at = 0;
            while(at < problem.T) {
                auto [l, v] = rec[at];
                dna[at].insert({v, l});
                int nxt = std::min(problem.T, at + problem.center_types[l]);
                facility_seg[v].update(at, nxt-1, -1);
                for (int t = at; t < nxt; t++) {
                    vis[t][v] = true;
                }
                at = std::min(problem.T, at + problem.center_types[l]);
            }
        }
        get_fitness();
    }


    int get_fitness() {
        fitness = 0;
        for (int t = 0; t < problem.T; t++) {
            for (auto [v, l]: dna[t]) {
                for (int i = t; i < std::min(problem.T, t + problem.center_types[l]); i++) {
                    current_facilities[i].insert(v);
                }
            }
        }

        for (int t = 0; t < problem.T; t++) {
            for (int client: problem.clients[t]) {
                long long mn = INF;
                for (int facility: current_facilities[t]) {
                    if (facility == -1) break;
                    mn = std::min(mn, problem.graph[client][facility]);
                }
                if (problem.version == "LKM") {
                    fitness += mn;
                }
                else if (problem.version == "LKC") {
                    fitness = std::max(fitness, mn);
                }
            }
        }
        
        return fitness;
    }
};