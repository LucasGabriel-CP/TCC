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
    int n_nodes, n_steiner, sub_connections;
    std::vector<int> go, energy, out_cable;
    float fitness=-1;

    Individuo() { }

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
        assert(con_sub > 0 && con_sub <= C);

        std::queue<int> que;
        for (int i = 0; i < n_nodes; i++) {
            if (!connections[i]) que.push(i);
        }
        assert(!que.empty());

        while (!que.empty()) {
            int u = que.front(); que.pop();
            out_cable[u] = best_cable(energy[u]).id;
            int v = go[u];
            if (v == NIL) continue;
            connections[v]--;
            energy[v] += energy[u];
            if (!connections[v]) que.push(v);
        }
        
        for (int &i: connections) assert(!i);
    }

    std::vector<std::vector<int>> create_layers(std::vector<int> const& order, int start) {
        std::vector<int> layer;
        std::vector<std::vector<int>> layers;
        int total_energy = 0;
        int biggest_capacity = cables.back().capacity;
        debug(start);
        for (int i = start - n_nodes; ;i++) {
            int at = order[i < 0 ? i + n_nodes : i];
            debug(i, at);
            if (i == start || total_energy + turbinies[at].total_prod > biggest_capacity) {
                debug(layer);
                layers.push_back(std::move(layer));
                total_energy = 0;
            }
            if (i == start) return layers;
            total_energy += turbinies[at].total_prod;
            layer.push_back(at);
        }
        assert(false);
    }

    void build(std::vector<int> const& order, int start) {
        std::vector<std::vector<int>> layers = create_layers(order, start);
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
        assert(validate());
        get_fitness();
    }

    void fix_this() {

    }

    void mutate() {
        float min_cost = FLOAT_INF;
        int min_i = -1, min_j = -1;
        for (int i = 0; i < n_nodes - 2; i++) {
            for (int j = i + 1; j < n_nodes; j++) {
                float change = (dist(i, go[i]) + dist(j, go[j])) - (dist(i, go[j]) + dist(j, go[i]));
                if (change < min_cost) {
                    min_i = i;
                    min_j = j;
                }
            }
        }

        if (f_cmp(min_cost, 0.0F)) {
            std::swap(go[min_i], go[min_j]);
            if (!validate()) {
                fix_this();
            }
            get_fitness();
        }
    }

    static std::pair<Individuo, Individuo> cross(Individuo const &p1, Individuo const &p2) {
        Individuo children_1(N), children_2(N);
        if (rand() < 0.5) {
            for (int i = 0; i < N/2; i++) {
                children_1.go[i] = p1.go[i];
                children_2.go[i] = p2.go[i];
                
                children_1.go[N-i-1] = p1.go[N-i-1];
                children_2.go[N-i-1] = p2.go[N-i-1];
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

        if (!children_1.validate()) {
            children_1.fix_this();
        }
        if (!children_2.validate()) {
            children_2.fix_this();
        }
        children_1.get_fitness();
        children_2.get_fitness();

        return {children_1, children_2};
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

    float get_fitness() {
        fitness = 0.0;

        propagate_energy();
        for (int i = 0; i < n_nodes; i++) {
            fitness += dist(i, go[i]) * best_cable_price(energy[i]);
        }
        
        return fitness;
    }
    
    friend bool operator < (const Individuo& lhs, const Individuo& rhs){
        return lhs.fitness < rhs.fitness;
    }
};
