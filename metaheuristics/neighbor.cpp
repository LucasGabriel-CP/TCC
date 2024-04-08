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
                    // if (check_type(t, old_l, l)) {
                        available_n2.push_back({t, v, old_l, l});
                    // }
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

    void first_improve_ls() {
        std::vector<std::set<std::pair<int, int>>> save_dna = dna;
        int op = 3;
        if (!op) {
            std::vector<ExchangeCenterType> available_n2;
            for (int t = 0; t < problem.T; t++) {
                for (auto [v, old_l]: dna[t]) {
                    for (int l = 0; l < problem.L; l++) {
                        if (l == old_l) continue;
                        available_n2.push_back({t, v, old_l, l});
                    }
                }
            }

            if (!available_n2.empty()) {
                int id = rand_i() % (int)available_n2.size();
                dna[available_n2[id].t].erase({available_n2[id].v, available_n2[id].l});
                dna[available_n2[id].t].insert({available_n2[id].v, available_n2[id].oth});
            }
        }
        else if (op == 1) {
            int t = rand_i() % problem.T, v = rand_i() % problem.V, l = rand_i() % problem.L;
            dna[t].insert({v, l});
        }
        else if (op == 2) {
            std::vector<ExchangeCenterType> available;
            for (int t = 0; t < problem.T; t++) {
                for (auto [v, l]: dna[t]) available.push_back({t, v, l, v});
            }
            if (!available.empty()) {
                int id = rand_i() % (int)available.size();
                dna[available[id].t].erase({available[id].v, available[id].l});
            }
        }
        else {
            std::vector<ExchangeCenterType> available_n3;
            for (int t = 0; t < problem.T; t++) {
                for (auto [old_v, l]: dna[t]) {
                    for (int v = 0; v < problem.V; v++) {
                        available_n3.push_back({t, old_v, l, v});
                    }
                }
            }

            if (!available_n3.empty()) {
                int id = rand_i() % (int)available_n3.size();
                dna[available_n3[id].t].erase({available_n3[id].v, available_n3[id].l});
                dna[available_n3[id].t].insert({available_n3[id].oth, available_n3[id].l});
            }
        }
        // if (give_penalties()) fix_solution();
        get_fitness();
    }

    void best_improve_ls() {
        bool improve = true;
        while (improve) {
            improve = false;
            std::vector<std::set<std::pair<int, int>>> save_dna = dna;
            long long save_fitness = fitness;
            int op = rand_i() % 3;
            if (!op) {
                std::vector<ExchangeCenterType> available_n2;
                for (int t = 0; t < problem.T; t++) {
                    for (auto [v, old_l]: dna[t]) {
                        for (int l = 0; l < problem.L; l++) {
                            available_n2.push_back({t, v, old_l, l});
                        }
                    }
                }

                if (!available_n2.empty()) {
                    int id = rand_i() % (int)available_n2.size();
                    dna[available_n2[id].t].erase({available_n2[id].v, available_n2[id].l});
                    dna[available_n2[id].t].insert({available_n2[id].v, available_n2[id].oth});
                }
            }
            else if (op == 1) {
                int t = rand_i() % problem.T, v = rand_i() % problem.V, l = rand_i() % problem.L;
                dna[t].insert({v, l});
            }
            else if (op == 2) {
                std::vector<ExchangeCenterType> available;
                for (int t = 0; t < problem.T; t++) {
                    for (auto [v, l]: dna[t]) available.push_back({t, v, l, v});
                }
                if (!available.empty()) {
                    int id = rand_i() % (int)available.size();
                    dna[available[id].t].erase({available[id].v, available[id].l});
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
            best_v[i] = {INF, i};
        }
        get_facilities();
        for (int t = 0; t < problem.T; t++) {
            for (int client: problem.clients[t]) {
                for (int i = 0; i < (int)solutions.size(); i++) {
                    auto sol = solutions[i];
                    int nd = std::min(problem.T, sol.t + problem.center_types[sol.l]);
                    if (sol.t > t || nd <= t) continue;
                    if (best_v[i].first == INF) best_v[i].first = 0;
                    best_v[i].first += problem.graph[client][sol.v];
                }
            }
        }

        std::sort(best_v.begin(), best_v.end());
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
            if ((int)current_facilities[t].empty()) need = true;
        }
        while(need) {
            int at = 0;
            for (at = 0; at < problem.T && !current_facilities[at].empty(); at++);
            int ct = -1, cv, cl;
            for (int t = at + 1; ct == -1 && t < problem.T; t++) {
                if (!dna[t].empty()) {
                    ct = t;
                    int op = rand_i() % (int)dna[t].size();
                    for (auto [v, l]: dna[t]) {
                        if (!op) {
                            cv = v; cl = l;
                            break;
                        }
                        op--;
                    }
                }
            }

            if (ct == -1) {
                cv = rand_i() % problem.V;
                cl = rand_i() % problem.L;
                dna[at].insert({cv, cl});
            }
            else {
                dna[at] = dna[ct];
                dna[ct].erase({cv, cl});
            }

            need = false;
            get_facilities();
            for (int t = 0; !need && t < problem.T; t++) {
                if ((int)current_facilities[t].empty()) need = true;
            }
        }
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