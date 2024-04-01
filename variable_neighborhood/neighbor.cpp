#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <assert.h>
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "util/string_hash.cpp"
#include "leasing.cpp"

struct Neighbor {
    struct ExchangeCenterType {
        int t, v, l, oth;
        ExchangeCenterType(int _t, int _v, int _l, int _oth) {
            t = _t;
            v = _v;
            l = _l;
            oth = _oth;
        }
    };

    long long fitness=-1;
    std::vector<std::set<std::pair<int, int>>> dna;
    std::vector<std::set<int>> current_facilities, facilities_types;
    std::vector<std::pair<int, int>> index_best;
    std::vector<long long> clients_best;
    LeasingProblem problem;

    Neighbor(LeasingProblem const &_problem) {
        problem = _problem;
        dna = std::vector<std::set<std::pair<int, int>>>(problem.T);
        current_facilities = std::vector<std::set<int>>(problem.T);
        facilities_types = std::vector<std::set<int>>(problem.T);
        index_best = std::vector<std::pair<int, int>>(problem.V);
        clients_best = std::vector<long long>(problem.V);
    }

    Neighbor(LeasingProblem const &_problem, std::vector<std::set<std::pair<int, int>>> _dna, long long _fitness) {
        problem = _problem;
        index_best = std::vector<std::pair<int, int>>(problem.V);
        clients_best = std::vector<long long>(problem.V);

        dna = _dna;
        current_facilities = std::vector<std::set<int>>(problem.T);
        facilities_types = std::vector<std::set<int>>(problem.T);
        fitness = _fitness;
        get_facilities();
    }

    bool check_facility(std::set<std::pair<int, int>> const &lhs, int t, int old_t) {
        for (auto [v, l]: lhs) {
            int nd = std::min(problem.T, t + problem.center_types[l]);
            int a_nd = std::min(problem.T, old_t + problem.center_types[l]);
            for (int st = t; st < nd; st++) {
                if (st >= old_t && st < old_t + a_nd) continue;
                if (current_facilities[st].find(v) != current_facilities[st].end()
                    || (int)lhs.size() + (int)current_facilities[st].size() > problem.K) {
                    return false;
                }
            }
            for (int lc: problem.center_types) {
                for (int st = std::max(0, t - lc + 1); st <= t; st++) {
                    if ((st >= old_t && st < old_t + a_nd) || facilities_types[st].find(lc) == facilities_types[st].end()) {
                        continue;
                    }
                    if (current_facilities[st].find(v) != current_facilities[st].end()) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool check_type(std::pair<int, int> const &lhs, int t, int l) {
        if (problem.center_types[l] <= problem.center_types[lhs.second]) {
            return true;
        }
        int nd = std::min(problem.T, t + problem.center_types[l]);
        for (int st = t + problem.center_types[lhs.second]; st < nd; st++) {
            if ((int)current_facilities[t].size() == problem.K) return false;
        }
        return true;
    }

    bool check_location(int t, int v) {
        for (int l: problem.center_types) {
            for (int st = std::max(0, t - l + 1); st < t; st++) {
                if (facilities_types[st].find(l) == facilities_types[st].end()) {
                    continue;
                }
                if (current_facilities[st].find(v) != current_facilities[st].end()) {
                    return false;
                }
            }
            for (int st = t; st < std::min(problem.T, t + l); st++) {
                if (facilities_types[st].find(l) == facilities_types[st].end()) {
                    continue;
                }
                if (current_facilities[st].find(v) != current_facilities[st].end()) {
                    return false;
                }
            }
        }
        return true;
    }

    bool shake() {
        bool shaked = false;
        get_facilities();
        std::set<std::pair<int, int>> available_n1;
        for (int t = 0; t < problem.T; t++) {
            for (int tt = t+1; tt < problem.T; tt++) {
                if (check_facility(dna[t], tt, t) && check_facility(dna[tt], t, tt)) {
                    available_n1.insert({t, tt});
                }
            }
        }
        if (!available_n1.empty()) {
            int id = rand_i() % (int)available_n1.size();
            int p1, p2;
            for (auto [t, tt]: available_n1) {
                if (!id) {
                    p1 = t; p2 = tt;
                    break;
                }
                --id;
            }
            std::swap(dna[p1], dna[p2]);
            get_fitness();
            shaked = true;
        }

        std::vector<ExchangeCenterType> available_n2;
        for (int t = 0; t < problem.T; t++) {
            for (auto [v, old_l]: dna[t]) {
                for (int l = 0; l < problem.L; l++) {
                    if (check_type({v, old_l}, t, l)) {
                        available_n2.push_back({t, v, old_l, l});
                    }
                }
            }
        }

        if (!available_n2.empty()) {
            int id = rand_i() % (int)available_n2.size();
            dna[available_n2[id].t].erase({available_n2[id].v, available_n2[id].l});
            dna[available_n2[id].t].insert({available_n2[id].v, available_n2[id].oth});
            get_fitness();
            shaked = true;
        }

        std::vector<ExchangeCenterType> available_n3;
        for (int t = 0; t < problem.T; t++) {
            for (auto [old_v, l]: dna[t]) {
                for (int v = 0; v < problem.V; v++) {
                    if (check_location(t, v)) {
                        available_n3.push_back({t, old_v, l, v});
                    }
                }
            }
        }

        if (!available_n3.empty()) {
            int id = rand_i() % (int)available_n3.size();
            dna[available_n3[id].t].erase({available_n3[id].v, available_n3[id].l});
            dna[available_n3[id].t].insert({available_n3[id].oth, available_n3[id].l});
            get_fitness();
            shaked = true;
        }
        return shaked;
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

        std::sort(index_best.begin(), index_best.end(), std::greater<std::pair<int, int>>());
        // std::shuffle(index_best.begin(), index_best.end(), g);
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
                i = nd-1;
            }
        }

        for (t = 0; t < problem.T; t++) {
            dna[t] = new_sol[t];
        }

        get_fitness();
    }

    void local_searchv2() {
        bool improve = true;
        while (improve) {
            improve = false;
            std::vector<std::set<std::pair<int, int>>> save_dna = dna;
            long long save_fitness = fitness;
            for (int t = 0; t < problem.T; t++) {
                int op = rand_i() % 3;
                if (!op) {
                    int v = rand_i() % problem.V;
                    int l = rand_i() % problem.L;
                    int nd = std::min(problem.T, t + problem.center_types[l]);
                    bool ok = true;
                    for (int st = t; ok && st < nd; st++) {
                        if ((int)current_facilities[st].size() == problem.K) {
                            ok = false;
                        }
                    }
                    if (check_location(t, v) && ok) {
                        dna[t].insert({v, l});
                        get_facilities();
                    }
                }
                else if (!dna[t].empty()) {
                    int rl, rv;
                    int k = rand_i() % (int)dna[t].size();
                    for (auto [i, j]: dna[t]) {
                        if (!k) {
                            rv = i; rl = j;
                            break;
                        }
                        k--;
                    }
                    int nd = std::min(problem.T, t + problem.center_types[rl]);
                    bool ok = true;
                    for (int st = t; ok && st < nd; st++) {
                        if ((int)current_facilities[st].size() == 1) {
                            ok = false;
                        }
                    }
                    if (ok) {
                        dna[t].erase({rv, rl});
                    }
                    int v = rand_i() % problem.V;
                    int l = rand_i() % problem.L;
                    if (op == 2) {
                        nd = std::min(problem.T, t + problem.center_types[l]);
                        ok = true;
                        for (int st = t; ok && st < nd; st++) {
                            if (current_facilities[st].find(v) != current_facilities[st].end()
                                || (int)current_facilities[st].size() == problem.K) {
                                ok = false;
                            }
                        }
                        if (check_location(t, v) && ok) {
                            dna[t].insert({v, l});
                            for (int st = t; st < nd; st++) {
                                current_facilities[st].insert(v);
                            }
                        }
                    }
                    get_facilities();
                }
            }
            if (get_fitness() < save_fitness) {
                improve = true;
            }
            else {
                dna = save_dna;
            }
        }
        get_fitness();
    }

    void get_facilities() {
        for (int t = 0; t < problem.T; t++) {
            current_facilities[t].clear();
            facilities_types[t].clear();
        }

        for (int t = 0; t < problem.T; t++) {
            for (auto [v, l]: dna[t]) {
                for (int i = t; i < std::min(problem.T, t + problem.center_types[l]); i++) {
                    current_facilities[i].insert(v);
                    facilities_types[i].insert(l);
                }
            }
        }
    }

    void valid() {
        for (int t = 0; t < problem.T; t++) {
            if ((int)current_facilities[t].size() > problem.K) {
                assert(false && "K constraint");
            }
            if (current_facilities[t].empty()) {
                assert(false && "No facility");
            }
            for (auto [v, l]: dna[t]) {
                int cnt = 0;
                int nd = std::min(problem.T, t + problem.center_types[l]);
                for (int st = t; st < nd; st++) {
                    // Check for multifacilities
                    for (auto [u, _]: dna[st]) {
                        if (u == v) {
                            cnt++;
                        }
                    }
                }
                if (cnt > 1) {
                    assert(false && "Multi facilityies");
                }
            }
        }
    }

    long long give_penalties() {
        long long total_penalty = 0, multi_facilities = 0, more_than_k = 0, no_facility = 0;

        for (int t = 0; t < problem.T; t++) {
            for (auto [v, l]: dna[t]) {
                int nd = std::min(problem.T, t + problem.center_types[l]);
                for (int st = t; st < nd; st++) {
                    // Check for multifacilities
                    for (auto [u, _]: dna[st]) {
                        if (u == v) {
                            multi_facilities++;
                        }
                    }
                }
                multi_facilities--;
            }
        }

        // Check for more than problem.K or no facility
        for (int t = 0; t < problem.T; t++) {
            if (!(int)current_facilities[t].size()) {
                no_facility++;
            }
            if ((int)current_facilities[t].size() > problem.K) {
                more_than_k++;
            }
        }

        total_penalty = multi_facilities * MULTI_FACILITY_PENALTY + more_than_k * PK_FACILITIES_PENALTY * PK_FACILITIES_PENALTY
                + no_facility * NO_FACILITY_PENALTY * NO_FACILITY_PENALTY * NO_FACILITY_PENALTY;

        return total_penalty;
    }

    long long get_fitness() {
        fitness = 0;

        get_facilities();
        // valid();
        fitness += give_penalties();
        for (int t = 0; t < problem.T; t++) {
            for (int client: problem.clients[t]) {
                long long mn = INF;
                for (int facility: current_facilities[t]) {
                    mn = std::min(mn, problem.graph[client][facility]);
                }
                if (mn == INF) continue;
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


    long long get_hashed_sol() {
        std::string sol="";
        for (int t = 0; t < problem.T; t++) {
            sol += ("T" + std::to_string(t));
            for (auto [v, l]: dna[t]) {
                sol += ("V" + std::to_string(v));
                sol += ("L" + std::to_string(l));
            }
        }
        HashedString hash(sol);
        return hash.get(0, (int)sol.size() - 1);
    }

    friend bool operator < (Neighbor const &lhs, Neighbor const &rhs){
        return lhs.fitness < rhs.fitness;
    }

    friend std::ostream &operator <<(std::ostream &os, const Neighbor &sol) {
        os << "Fitness: " << sol.fitness;
        return os;
    }
};