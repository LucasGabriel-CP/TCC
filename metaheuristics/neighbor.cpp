#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <assert.h>
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "util/string_hash.cpp"
#include "problems/leasing.cpp"

struct Neighbor {
    struct ExchangeCenterType {
        int t, v, l, oth;
        ExchangeCenterType() { }
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

    bool check_type(int t, int old_l, int l) {
        if (problem.center_types[l] <= problem.center_types[old_l]) {
            return true;
        }
        int nd = std::min(problem.T, t + problem.center_types[l]);
        for (int st = t + problem.center_types[old_l]; st < nd; st++) {
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

    bool shake(int r=10) {
        bool shaked = false;

        get_facilities();

        while (r--) {
            std::vector<ExchangeCenterType> available_n2;
            for (int t = 0; t < problem.T; t++) {
                for (auto [v, old_l]: dna[t]) {
                    for (int l = 0; l < problem.L; l++) {
                        if (l != old_l) {
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
        }
        
        return shaked;
    }

    void get_neighbor(int r = 3) {
        int cnt = 0;
        for (int t = 0; t < problem.T; t++) {
            cnt += (int)dna[t].size();
        }
        std::vector<std::pair<int, int>> best_v(cnt);
        std::vector<ExchangeCenterType> solutions(cnt), new_solution;
        for (int t = 0, i = 0; t < problem.T; t++) {
            for (auto [v, l]: dna[t]) {
                solutions[i++] = {t, v, l, -1};
            }
        }
        for (int i = 0; i < (int)solutions.size(); i++) {
            best_v[i] = {0, i};
        }
        get_facilities();
        for (int t = 0; t < problem.T; t++) {
            for (int client: problem.clients[t]) {
                long long mn = INF;
                int id = -1;
                for (int i = 0; i < (int)solutions.size(); i++) {
                    auto sol = solutions[i];
                    int nd = std::min(problem.T, sol.t + problem.center_types[sol.l]);
                    if (sol.t > t || nd <= t) continue;
                    if (problem.graph[client][sol.v] < mn) {
                        mn = problem.graph[client][sol.v];
                        id = i;
                    }
                }
                if (id != -1) {
                    best_v[id].first++;
                }
            }
        }

        // std::shuffle(best_v.begin(), best_v.end(), g);
        std::sort(best_v.begin(), best_v.end(), [&](auto a, auto b){ return a.first > b.first; });
        r = std::max(cnt - r, 1);
        for (int i = 0; i < r; i++) {
            new_solution.push_back(solutions[best_v[i].second]);
        }
        for (int t = 0; t < problem.T; t++) {
            dna[t].clear();
            current_facilities[t].clear();
            facilities_types[t].clear();
        }
        for (auto [t, v, l, _]: new_solution) {
            int nd = std::min(problem.T, t + problem.center_types[l]);
            for (int st = t; st < nd; st++) {
                current_facilities[t].insert(v);
                facilities_types[t].insert(l);
            }
        }

        bool filled = false;
        while(!filled){
            int t = 0, l = rand_i() % problem.L;
            int tt=t;
            filled = true;
            while(tt < problem.T){
                int st = tt;
                int nd = std::min(problem.T-1, st + problem.center_types[l] - 1);
                int sz = 0;
                for(int i = st; i <= nd; i++){
                    sz = std::max(sz, (int)current_facilities[i].size());
                }
                if(sz+1 <= problem.K){
                    std::vector<int> options_facilities;
                    for(int v = 0; v < problem.V; v++){
                        bool has = false;
                        for(int i = st; !has && i <= nd; i++){
                            if(current_facilities[i].find(v) != current_facilities[i].end()){
                                has=true;
                            }
                        }
                        if(!has){
                            options_facilities.push_back(v);
                        }
                    }
                    if(!options_facilities.empty()) {
                        int op = rand_i() % options_facilities.size();
                        int v = options_facilities[op];
                        for(int u = st; u <= nd; u++){
                            current_facilities[u].insert(v);
                        }
                        new_solution.push_back({tt, v, l, -1});
                        filled=false;
                        break;
                    }
                }
                tt++;
            }
        }
        for (int t = 0; t < problem.T; t++) {
            dna[t].clear();
        }
        for (auto sol: new_solution) {
            dna[sol.t].insert({sol.v, sol.l});
        }
        get_fitness();
    }

    void first_improve_ls(int r=3) {
        std::vector<std::set<std::pair<int, int>>> save_dna = dna;
        long long save_fitness = fitness;
        int cnt = 0;
        for (int t = 0; t < problem.T; t++) {
            cnt += (int)dna[t].size();
        }
        do {
            save_dna = dna;
            std::vector<ExchangeCenterType> available_n3;
            save_fitness = fitness;
            for (int t = 0; t < problem.T; t++) {
                for (auto [old_v, l]: dna[t]) {
                    for (int v = 0; v < problem.V; v++) {
                        if (v != old_v) {
                            available_n3.push_back({t, old_v, l, v});
                        }
                    }
                }
            }
            std::shuffle(available_n3.begin(), available_n3.end(), g);

            long long mn_fitness = fitness;
            std::vector<std::set<std::pair<int, int>>> best_dna = dna;
            int nd = std::min(cnt, (int)available_n3.size());
            for (int id = 0; id < nd; id++) {
                dna[available_n3[id].t].erase({available_n3[id].v, available_n3[id].l});
                dna[available_n3[id].t].insert({available_n3[id].oth, available_n3[id].l});
                if (give_penalties()) fix_solution();
                get_fitness();
                if (fitness > mn_fitness) {
                    dna = best_dna;
                    fitness = mn_fitness;
                }
                else {
                    best_dna = dna;
                    mn_fitness = fitness;
                }
            }
        } while(fitness >= save_fitness && r--);
    }

    void best_improve_ls() {
        bool improve = true;
        while (improve) {
            improve = false;
            std::vector<std::set<std::pair<int, int>>> save_dna = dna;
            long long save_fitness = fitness;
            std::vector<ExchangeCenterType> available_n3;
            for (int t = 0; t < problem.T; t++) {
                for (auto [old_v, l]: dna[t]) {
                    for (int v = 0; v < problem.V; v++) {
                        if (v != old_v) {
                            available_n3.push_back({t, old_v, l, v});
                        }
                    }
                }
            }
            std::shuffle(available_n3.begin(), available_n3.end(), g);

            long long mn_fitness = fitness;
            std::vector<std::set<std::pair<int, int>>> best_dna = dna;
            for (int id = 0; id < (int)available_n3.size(); id++) {
                dna[available_n3[id].t].erase({available_n3[id].v, available_n3[id].l});
                dna[available_n3[id].t].insert({available_n3[id].oth, available_n3[id].l});
                if (give_penalties()) fix_solution();
                get_fitness();
                if (fitness > mn_fitness) {
                    dna = best_dna;
                    fitness = mn_fitness;
                }
                else {
                    // available_n3.clear();
                    // for (int t = 0; t < problem.T; t++) {
                    //     for (auto [old_v, l]: dna[t]) {
                    //         for (int v = 0; v < problem.V; v++) {
                    //             if (v != old_v) {
                    //                 available_n3.push_back({t, old_v, l, v});
                    //             }
                    //         }
                    //     }
                    // }
                    // std::shuffle(available_n3.begin(), available_n3.end(), g);
                    best_dna = dna;
                    mn_fitness = fitness;
                }
            }
            if (get_fitness() < save_fitness) {
                improve = true;
            }
            else {
                dna = save_dna;
                fitness = save_fitness;
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
        #ifdef DPIZZA
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
        #endif
    }

    long long give_penalties() {
        long long total_penalty = 0, multi_facilities = 0, more_than_k = 0, no_facility = 0;

        std::set<std::pair<int, int>> used;

        for (int t = 0; t < problem.T; t++) {
            for (auto [v, l]: dna[t]) {
                int nd = std::min(problem.T, t + problem.center_types[l]);
                int added = 0;
                for (int st = t; st < nd; st++) {
                    // Check for multifacilities
                    if (used.find({st, v}) != used.end()) {
                        added++;
                    }
                    used.insert({st, v});
                }
                if (added)
                    multi_facilities += added;
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
        // long long penalties = give_penalties();
        // if (penalties != 0) {
        //     fix_solution();
        // }
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

        fitness += give_penalties();
        
        return fitness;
    }

    void fix_solution() {
        std::vector<std::set<std::pair<int, int>>> fixed_sol(problem.T);
        int cnt = 0;
        for (int t = 0; t < problem.T; t++) {
            cnt += (int)dna[t].size();
        }
        std::vector<std::pair<int, int>> best_v(cnt);
        std::vector<ExchangeCenterType> solutions(cnt);
        for (int t = 0, i = 0; t < problem.T; t++) {
            for (auto [v, l]: dna[t]) {
                solutions[i++] = {t, v, l, -1};
            }
        }
        for (int i = 0; i < (int)solutions.size(); i++) {
            best_v[i] = {0, i};
        }
        get_facilities();
        for (int t = 0; t < problem.T; t++) {
            for (int client: problem.clients[t]) {
                long long mn = INF;
                int id = -1;
                for (int i = 0; i < (int)solutions.size(); i++) {
                    auto sol = solutions[i];
                    int nd = std::min(problem.T, sol.t + problem.center_types[sol.l]);
                    if (sol.t > t || nd <= t) continue;
                    if (problem.graph[client][sol.v] < mn) {
                        mn = problem.graph[client][sol.v];
                        id = i;
                    }
                }
                if (id != -1) {
                    best_v[id].first++;
                }
            }
        }

        // std::shuffle(best_v.begin(), best_v.end(), g);
        std::sort(best_v.begin(), best_v.end(), [&](auto a, auto b){ return a.first > b.first; });
        for (int t = 0; t < problem.T; t++) {
            dna[t].clear();
            current_facilities[t].clear();
            facilities_types[t].clear();
        }
        for (auto [_, i]: best_v) {
            auto sol = solutions[i];
            bool ok = true;
            for (int st = sol.t; ok && st < std::min(problem.T, sol.t + problem.center_types[sol.l]); st++) {
                if ((int)current_facilities[st].size () == problem.K || current_facilities[st].find(sol.v) != current_facilities[st].end()) {
                    ok = false;
                }
            }
            if (ok) {
                dna[sol.t].insert({sol.v, sol.l});
                for (int st = sol.t; st < std::min(problem.T, sol.t + problem.center_types[sol.l]); st++) {
                    current_facilities[st].insert(sol.v);
                    facilities_types[st].insert(sol.l);
                }
            }
        }
        bool need = false;
        for (int t = 0; !need && t < problem.T; t++) {
            if ((int)current_facilities[t].size() < problem.K) need = true;
            assert((int)current_facilities[t].size() <= problem.K);
        }
        int at = 0;
        while(true) {
            for (; at < problem.T && (int)current_facilities[at].size() == problem.K; at++);
            if (at >= problem.T) break;

            std::set<std::pair<int, int>> available;
            for (int v = 0; v < problem.V; v++) {
                for (int l = 0; l < problem.L; l++) {
                    int nd = std::min(problem.T, at + problem.center_types[l]), st;
                    for (st = at; st < nd && (int)current_facilities[st].size() < problem.K && current_facilities[st].find(v) == current_facilities[st].end(); st++);
                    if (st == nd) {
                        available.insert({v, l});
                    }
                }
            }
            if (available.empty()) {
                if (current_facilities[at].empty()) {
                    int st;
                    for (st = at + 1; st < problem.T && dna[st].empty(); st++);
                    assert(st < problem.T);
                    std::swap(dna[at], dna[st]);
                    get_facilities();
                }
                at++;
            }
            else {
                int op = rand_i() % (int)available.size();
                int cv, cl;
                for (auto [v, l]: available) {
                    if (!op) {
                        cv = v; cl = l; break;
                    }
                    op--;
                }

                dna[at].insert({cv, cl});

                get_facilities();
            }
        }
        valid();
    }


    long long get_hashed_sol() {
        std::string sol="";
        for (int t = 0; t < problem.T; t++) {
            if (!dna[t].empty())
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