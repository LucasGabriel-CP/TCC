#pragma once

#include <iostream>
#include <random>
#include <algorithm>
#include <set>
#include <assert.h>

#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "problems/leasing.cpp"


struct Individuo {
    long long fitness=-1;
    std::vector<std::set<std::pair<int, int>>> dna;
    std::vector<std::set<int>> current_facilities, facilities_types;
    std::vector<long long> clients_best;
    std::vector<std::pair<int, int>> index_best;
    LeasingProblem problem;

    struct ExchangeCenterType {
        int t, v, l, oth;
        ExchangeCenterType(int _t, int _v, int _l, int _oth) {
            t = _t;
            v = _v;
            l = _l;
            oth = _oth;
        }
    };

    Individuo() {
        dna = std::vector<std::set<std::pair<int, int>>>(problem.T);
        clients_best = std::vector<long long>(problem.V, INF);
        index_best = std::vector<std::pair<int, int>>(problem.T*problem.K);
        current_facilities = std::vector<std::set<int>>(problem.T);
        facilities_types = std::vector<std::set<int>>(problem.T);
    }

    Individuo(LeasingProblem const &_problem) {
        problem = _problem;
        dna = std::vector<std::set<std::pair<int, int>>>(problem.T);
        clients_best = std::vector<long long>(problem.V, INF);
        index_best = std::vector<std::pair<int, int>>(problem.T*problem.K);
        current_facilities = std::vector<std::set<int>>(problem.T);
        facilities_types = std::vector<std::set<int>>(problem.T);
    }

    Individuo(std::vector<std::set<std::pair<int, int>>> _matrix, int _fitness, LeasingProblem const &_problem) {
        problem = _problem;
        dna = _matrix;
        fitness = _fitness;
        current_facilities = std::vector<std::set<int>>(problem.T);
        facilities_types = std::vector<std::set<int>>(problem.T);
        clients_best = std::vector<long long>(problem.V, INF);
        index_best = std::vector<std::pair<int, int>>(problem.T*problem.K);
        get_facilities();
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

    void get_facility_LKM(int t) {
        int bv = -1, bl = -1, orig = t;
        long long best = INF;
        while (bl == -1) {
            for (int l = 0; l < problem.L; l++) {
                int nd = std::min(problem.T, t + problem.center_types[l]);
                assert(nd > orig);
                bool ok = true;
                for (int st = t; st < nd; st++) {
                    if ((int)current_facilities[st].size() == problem.K) ok = false;
                }
                if (!ok) continue;
                for (int v = 0; v < problem.V; v++) {
                    ok = true;
                    for (int st = t; st < nd; st++) {
                        if (current_facilities[st].find(v) != current_facilities[st].end()) {
                            ok = false;
                        }
                    }
                    if (!ok) continue;
                    long long sum = 0;
                    for (int client: problem.clients[t]) {
                        sum += problem.graph[client][v];
                    }
                    if (best > sum) {
                        bv = v; bl = l;
                        best = sum;
                    }
                }
            }
            if (bl == -1) --t;
        }
        dna[t].insert({bv, bl});
        get_facilities();
    }

    // void mutate() {
    //     // IMPROVE THIS LIKE THE SHAKE ON VNS
    //     for (int t = 0; t < problem.T; t++) {
    //         current_facilities[t].clear();
    //         dna[t].clear();
    //     }
    //     for (int t = 0; t < problem.T; t++) {
    //         // if ((int)current_facilities[t].size() == problem.K) continue;
    //         if (rand() < problem.K / (problem.mean_l * problem.V)) {
    //             for (int v = 0; v < problem.V; v++) {
    //                 int l = rand_i() % problem.L;
    //                 int nd = std::min(problem.T, t + problem.center_types[l]);
    //                 // bool ok = true;
    //                 for (int st = t; st < nd; st++) {
    //                     current_facilities[st].insert(v);
    //                     facilities_types[st].insert(l);
    //                 }
    //                 dna[t].insert({v, l});
    //             }
    //         }
    //     }

    //     get_fitness();
    // }

    void mutate() {
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
        }

        std::vector<ExchangeCenterType> available_n2;
        for (int t = 0; t < problem.T; t++) {
            for (auto [v, old_l]: dna[t]) {
                for (int l = 0; l < problem.L; l++) {
                    // if (check_type({v, old_l}, t, l)) {
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
        }

        std::vector<ExchangeCenterType> available_n3;
        for (int t = 0; t < problem.T; t++) {
            for (auto [old_v, l]: dna[t]) {
                for (int v = 0; v < problem.V; v++) {
                    // if (check_location(t, v)) {
                        available_n3.push_back({t, old_v, l, v});
                    // }
                }
            }
        }

        if (!available_n3.empty()) {
            int id = rand_i() % (int)available_n3.size();
            dna[available_n3[id].t].erase({available_n3[id].v, available_n3[id].l});
            dna[available_n3[id].t].insert({available_n3[id].oth, available_n3[id].l});
            get_fitness();
        }
    }

    bool try_insert(std::set<std::pair<int, int>> gene, Individuo &child, int t) {
        for (auto [v, l]: gene) {
            int nd = std::min(problem.T, problem.center_types[l] + t);
            child.dna[t].insert({v, l});
            for (int st = t; st < nd; st++) {
                child.current_facilities[st].insert(v);
                child.facilities_types[st].insert(l);
            }
        }
        return true;
    }

    bool try_insert(std::pair<int, int> alele, Individuo &child, int t) {
        auto [v, l] = alele;
        int nd = std::min(problem.T, problem.center_types[l] + t);
        child.dna[t].insert({v, l});
        for (int st = t; st < nd; st++) {
            child.current_facilities[st].insert(v);
            child.facilities_types[st].insert(l);
        }
        return true;
    }

    void cross(Individuo const &p1, Individuo const &p2) {
        Individuo children_1(problem), children_2(problem);
        if (rand() < 0.5) {
            int cnt = 0;
            for (int t = 0; t < problem.T; t++) {
                cnt += std::max((int)p1.dna[t].size(), (int)p2.dna[t].size());
            }
            int lim = rand_i() % std::max(1, cnt);
            for (int t = 0, i = 0; t < problem.T; t++) {
                bool ok_1 = false, ok_2 = false;
                for (auto alele: p1.dna[t]) {
                    if (i < lim) {
                        ok_1 |= try_insert(alele, children_1, t);
                    }
                    else {
                        ok_2 |= try_insert(alele, children_2, t);
                    }
                    i++;
                }
                for (auto alele: p2.dna[t]) {
                    if (i < lim) {
                        ok_2 |= try_insert(alele, children_2, t);
                    }
                    else {
                        ok_1 |= try_insert(alele, children_1, t);
                    }
                    i++;
                }
            }
        }
        else {
            int l = rand_i() % problem.T, r = rand_i() % problem.T;
            for (int t = 0; t < problem.T; t++) {
                if (t >= l && t <= r) {
                    try_insert(p2.dna[t], children_1, t);
                    try_insert(p1.dna[t], children_2, t);
                }
                else {
                    try_insert(p1.dna[t], children_1, t);
                    try_insert(p2.dna[t], children_2, t);
                }
            }
        }

        // Select a random child to return
        if (children_1.get_fitness() < children_2.get_fitness()) {
            dna = children_1.dna;
        }
        else {
            dna = children_2.dna;
        }
        get_fitness();
    }

    bool validate() {        
        #ifdef PIZZA
            // Verify solution
        #endif

        return true;
    }


    void local_search() {
        Individuo new_sol(problem);
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
        int choosen = rand_i() % cnt; 

        for (int i = 0; i < choosen; i++) {
            int id = index_best[i].second;
            for (int t = 0, c = 0; t < problem.T; t++) {
                if (dna[t].empty()) continue;
                if (c + (int)dna.size() <= id) {
                    for (auto [u, l]: dna[t]) {
                        if (c == id) {
                            new_sol.dna[t].insert({u, l});
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
                    new_sol.dna[i].insert({v, l});
                    ok = true;
                    t = i+1;
                    break;
                }
                i = nd;
            }
        }

        for (t = 0; t < problem.T; t++) {
            dna[t] = new_sol.dna[t];
        }

        get_fitness();
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
                    // int nd = std::min(problem.T, t + problem.center_types[l]);
                    // bool ok = true;
                        dna[t].insert({v, l});
                        get_facilities();
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
                    // bool ok = true;
                        dna[t].erase({rv, rl});
                    int v = rand_i() % problem.V;
                    int l = rand_i() % problem.L;
                    if (op == 2) {
                        nd = std::min(problem.T, t + problem.center_types[l]);
                        // ok = true;
                            dna[t].insert({v, l});
                            for (int st = t; st < nd; st++) {
                                current_facilities[st].insert(v);
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

        for (int t = 0; t < problem.T; t++) {
            for (int client: problem.clients[t]) {
                long long mn = INF;
                for (int facility: current_facilities[t]) {
                    mn = std::min(mn, problem.graph[client][facility]);
                }
                if (mn != INF) {
                    if (problem.version == "LKC") {
                        fitness = std::max(fitness, mn);
                    }
                    else if (problem.version == "LKM") {
                        fitness += mn;
                    }
                }
            }
        }
        
        fitness += give_penalties();
        
        return fitness;
    }
    
    friend bool operator < (const Individuo &lhs, const Individuo &rhs){
        return lhs.fitness < rhs.fitness;
    }

    friend std::ostream &operator <<(std::ostream &os, const Individuo &ind) {
        os << "Fitness: " << ind.fitness;
        return os;
    }
};