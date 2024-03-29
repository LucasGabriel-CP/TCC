#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <queue>

#include "util/Cable.cpp"
#include "util/Point.cpp"
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "util/Turbine.cpp"
#include "SWFCR.cpp"

struct Individuo {
    int n_nodes, n_steiner, sub_connections, sub_energy;
    std::vector<int> go, energy, out_cable;
    double fitness=-1;

    Individuo() {
        n_nodes = n_steiner = sub_connections = sub_energy = -1;
        go = energy = out_cable = {};
    }

    Individuo(int _n_nodes) {
        n_nodes = _n_nodes;
        n_steiner = 0;
        go = std::vector<int>(n_nodes);
        out_cable = std::vector<int>(n_nodes);
    }

    Individuo(int _n_nodes, int _n_steiner) {
        n_nodes = _n_nodes;
        n_steiner = _n_steiner;
        out_cable = std::vector<int>(n_nodes);
        go = std::vector<int>(n_nodes);
    }

    void propagate_energy() {
        std::vector<int> connections(n_nodes, 0);
        energy = std::vector<int>(n_nodes, 1);
        int con_sub = 0;
        for (int i = 0; i < n_nodes; i++) {
            if (go[i] != NIL) {
                connections[go[i]]++;
            }
            else {
                con_sub++;
            }
        }
        // assert(con_sub > 0 && con_sub <= C);

        std::queue<int> que;
        for (int i = 0; i < n_nodes; i++) {
            if (!connections[i]) que.push(i);
        }
        // assert(!que.empty());

        while (!que.empty()) {
            int u = que.front(); que.pop();
            out_cable[u] = best_cable(energy[u]).id;
            int v = go[u];
            if (v == NIL) {
                sub_energy += energy[u];
            }
            else {
                connections[v]--;
                energy[v] += energy[u];
                if (!connections[v]) que.push(v);
            }
        }
        
        // for (int &i: connections) assert(!i);
    }

    std::vector<std::vector<int>> create_layers(std::vector<int> const& order, int start, bool ccw = false) {
        std::vector<int> layer;
        std::vector<std::vector<int>> layers;
        int total_energy = 0;
        int biggest_capacity = cables.back().capacity;
        bool setted = false;
        for (int i = start; ; (ccw ? i-- : i++)) {
            int at = order[(i + n_nodes) % n_nodes];
            if ((i + n_nodes) % n_nodes == start) setted = !setted;
            if (!setted || total_energy + turbinies[at].total_prod > biggest_capacity) {
                layers.push_back(std::move(layer));
                total_energy = 0;
            }
            if (!setted) return layers;
            total_energy += turbinies[at].total_prod;
            layer.push_back(at);
        }
        assert(false);
    }

    void build(std::vector<int> const& order, int start, bool ccw = false) {
        std::vector<std::vector<int>> layers = create_layers(order, start, ccw);
        if (ccw) {
            mtx.lock();
            debug(layers);
            mtx.unlock();
        }
        int cnt = 0;
        for (std::vector<int> layer: layers) {
            int u = layer.back();
            cnt += layer.size();
            for (int i = (int)layer.size() - 2; i >= 0; i--) {
                go[u] = layer[i];
                u = layer[i];
            }
            go[u] = NIL;
            sub_connections++;
        }
        assert(cnt == N);
        // assert(validate());
        get_fitness();
    }

    void fix_this() {

    }

    void mutate() {
        double min_cost = fitness;
        int min_i = -1, min_j = -1;
        for (int i = 0; i < n_nodes - 1; i++) {
            int aux = go[i];
            for (int j = i + 1; j < n_nodes; j++) {
                go[i] = j;
                double change = get_fitness();
                if (change < min_cost) {
                    min_i = i;
                    min_j = j;
                    min_cost = change;
                }
                go[i] = aux;
            }
        }

        if (min_i != -1) {
            go[min_i] = min_j;
        }
        else {
            int u = rand_i() % n_nodes;
            int v = rand_i() % (n_nodes + 1);
            while (v == u) v = rand_i() % (n_nodes + 1);
            if (v == n_nodes) v = NIL;
            go[u] = v;
        }

        get_fitness();
    }

    static Individuo cross(Individuo const &p1, Individuo const &p2) {
        return p1;
        Individuo children_1(N), children_2(N);
        assert((int)p1.go.size() == N);
        assert((int)p2.go.size() == N);
        if (rand() < 0.5) {
            for (int i = 0; i < N/2; i++) {
                children_1.go[i] = p1.go[i];
                children_2.go[i] = p2.go[i];
                
                children_2.go[N-i-1] = p1.go[N-i-1];
                children_1.go[N-i-1] = p2.go[N-i-1];
            }
        }
        else {
            for (int i = 0; i < N; i++) {
                if (rand() < 0.5) {
                    children_1.go[i] = p1.go[i];
                    children_2.go[i] = p2.go[i];
                }
                else {
                    children_2.go[i] = p1.go[i];
                    children_1.go[i] = p2.go[i];                                        
                }
            }
        }

        // if (!children_1.validate()) {
        //     children_1.fix_this();
        // }
        // if (!children_2.validate()) {
        //     children_2.fix_this();
        // }
        children_1.get_fitness();
        children_2.get_fitness();
        if (rand() < 0.5) {
            return children_1;
        }
        return children_2;
    }

    bool validate() {        
        #ifdef PIZZA
            for (int i = 0; i < n_nodes; i++) {
                assert(go[i] != i);
            }
        #endif

        int sub_c = 0;
        std::vector<char> vis(n_nodes, 'w');
        for (int u = 0; u < n_nodes; u++) {
            if (go[u] == NIL) sub_c++;
            if (vis[u] == 'w') {
                int v;
                for (v = u; v != NIL && vis[v] == 'w'; v = go[v]) {
                    vis[v] = 'g';
                }
                if (v != NIL && vis[v] == 'g') return false;
                for (v = u; v != NIL && vis[v] == 'g'; v = go[v]) {
                    vis[v] = 'b';
                }
            }
        }

        if (!sub_c || sub_c > C) return false;

        return true;
    }

    double give_penalties() {
        int cicle_penalties = 0, sub_c = 0, disconnected = 0, crossed = 0;
        bool not_reach = 0, C_overflow = 0;
        std::vector<char> vis(n_nodes, 'w');
        for (int u = 0; u < n_nodes; u++) {
            if (go[u] == NIL) sub_c++;
            // if (go[go[u]] == u) cicle_penalties++;
            if (vis[u] == 'w') {
                int v;
                for (v = u; v != NIL && vis[v] == 'w'; v = go[v]) {
                    vis[v] = 'g';
                }
                if (v != NIL && vis[v] == 'g') return cicle_penalties++;
                for (v = u; v != NIL && vis[v] == 'g'; v = go[v]) {
                    vis[v] = 'b';
                }
            }
        }
        not_reach = !sub_c;
        C_overflow = sub_c > C;

        std::vector<bool> dp(n_nodes, false);
        vis = std::vector<char>(n_nodes, 'w');
        std::function<bool(int)> dfs = [&](int u) {
            if (u == NIL) return true;
            if (vis[u] != 'w') return (bool)dp[u];
            vis[u] = 'g';
            return (bool)(dp[u] = dfs(go[u]));
        };

        for (int u = 0; u < n_nodes; u++) {
            if (!dfs(u)) disconnected++;
        }



        for (int u = 0; u < n_nodes-1; u++) {
            Turbine p = turbinies[u], q = (go[u] == NIL ? Turbine(substation) : turbinies[go[u]]);
            for (int v = u+1; v < n_nodes; v++) {
                if (go[u] == v || go[v] == u) continue;
                Turbine r = turbinies[v], s = (go[v] == NIL ? Turbine(substation) : turbinies[go[v]]);
                crossed += intersect(p, q, r, s);
            }
        }

        double total_penalty = 0.0;
        total_penalty += cicle_penalties * CLICE_PENALTY;
        total_penalty += disconnected * DISCONNECTED_PENALTY;
        total_penalty += not_reach * NOT_REACH_PENALTY;
        total_penalty += C_overflow * C_OVERFLOW_PENALTY;
        total_penalty += crossed * CROSS_PENLATY;

        return total_penalty;
    }

    double get_fitness() {
        fitness = 0.0;

        propagate_energy();
        for (int i = 0; i < n_nodes; i++) {
            fitness += dist(i, go[i]) * best_cable_price(energy[i]);
        }

        fitness += give_penalties();
        
        return fitness;
    }
    
    friend bool operator < (const Individuo lhs, const Individuo rhs){
        return lhs.fitness < rhs.fitness;
    }

    friend std::ostream &operator <<(std::ostream &os, const Individuo &ind) {
        os << "Fitness: " << ind.fitness << '\n' << "Energy: ";
        for (int i = 0; i < ind.n_nodes; i++) {
            os << ind.energy[i] << " ";
        }
        os << "\nSending: ";
        for (int i = 0; i < ind.n_nodes; i++) {
            os << ind.go[i] << " ";
        }
        os << "\nCables: ";
        for (int i = 0; i < ind.n_nodes; i++) {
            os << ind.out_cable[i] << " ";
        }
        return os;
    }
};
