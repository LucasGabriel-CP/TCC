#include <iostream>
#include <thread>
#include <assert.h>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
#include <set>
#include <unordered_map>
#include "util/Constants.cpp"
#include "util/DebugTemplate.cpp"

int V, T, L, K, E;
std::vector<std::vector<std::pair<int, int>>> graph, ordered_graph;      // dist(i, j)
std::vector<std::vector<int>> clients;          // Clientes no instante t
std::vector<int> center_types;                  // Tipos de centro

void read(std::string const &instance_path) {
    std::ifstream in;

    debug(instance_path);
    in.open(instance_path);

    in >> V >> T >> L >> K;
    debug(V, T, L, K);

    graph.assign(V, std::vector<std::pair<int, int>>(V));
    ordered_graph.assign(V, std::vector<std::pair<int, int>>(V));
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            graph[i][j].first = j;
            in >> graph[i][j].second;
            ordered_graph[i][j] = graph[i][j];
        }
        std::sort(ordered_graph[i].begin(), ordered_graph[i].end(), [](auto a, auto b){
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
    }
    debug(clients);

    in.close();
}



struct SegTree{
    #define Mid ((l + r) >> 1)
    #define nil 0
    #define MAXT 155
    int n;
    int lazy[4*MAXT], seg[4*MAXT];
    int modify(int &a, int &b){
        return int(a + b);
    }
    void push(int id, int tam){
        if (lazy[id]){
            seg[id] += lazy[id] * (tam + 1);
            if (tam){
                lazy[2*id] += lazy[id];
                lazy[2*id+1] += lazy[id];
            }
            lazy[id] = 0;
        }
    }
    int query(int id, int l, int r, int x, int y){
        push(id, r-l);
        if(x > r || y < l)   return nil; //doens'int change ans
        if(l >= x && r <= y) return seg[id];
        int p1 = query(2*id, l, Mid, x, y);
        int p2 = query(2*id+1, Mid+1, r, x, y);
        return modify(p1, p2);
    }
    void update(int id, int l, int r, int x, int y, int val){
        push(id, r-l);
        if(x > r || y < l)  return;
        if(l >= x && r <= y){
            lazy[id] = val;
            push(id, r-l);
            return;
        }
        update(2*id, l, Mid, x, y, val);
        update(2*id+1, Mid+1, r, x, y, val);
        seg[id] = modify(seg[2*id], seg[2*id+1]);
    }
    SegTree() {
        n = T;
        for (int i = 0; i < 4 * n; i++) {
            lazy[i] = 0;
            seg[i] = 0;
        }
    }
    SegTree(int _n){
        n = _n;
        for (int i = 0; i < 4 * n; i++) {
            lazy[i] = 0;
            seg[i] = 0;
        }
    }
    int query(int x, int y){ return query(1, 0, n-1, x, y); }
    void update(int x, int y, int val){ update(1, 0, n-1, x, y, val); }
};


struct InitialSolution {
    int fitness=-1;
    std::vector<std::set<std::pair<int, int>>> dna;
    std::vector<std::vector<int>> current_facilities;
    std::vector<std::pair<int, int>> rec;
    std::vector<std::vector<bool>> vis;
    std::vector<int> aproximation_server, dp, open_facilities, facility_idx;
    std::vector<SegTree> facility_seg;


    InitialSolution() {
        dna = std::vector<std::set<std::pair<int, int>>>(T);
        current_facilities = std::vector<std::vector<int>>(T, std::vector<int>(V, -1));
        rec = std::vector<std::pair<int, int>>(T);
        aproximation_server = std::vector<int>(K, -1);
        dp = std::vector<int>(T, -1);
        open_facilities = std::vector<int>(T, 0);
        facility_idx = std::vector<int>(T, 0);
        vis = std::vector<std::vector<bool>>(T, std::vector<bool>(V, false));
        facility_seg = std::vector<SegTree>(V, SegTree(T));
    }

    
    void dorit(int t) {
        int low = 0, high = V-1;
        std::vector<int> ans;
        while (low < high) {
            int mid = (low + high) / 2;
            std::vector<int> S;
            std::set<int> set_V;
            for (int i = 0; i < V; i++) set_V.insert(i);
            while (!set_V.empty()) {
                int x = *set_V.begin();
                S.push_back(x);
                for (int i = 0; i <= mid; i++) {
                    int v = ordered_graph[x][i].first;
                    if (set_V.find(v) != set_V.end()) {
                        set_V.erase(v);
                        for (int j = 0; j <= mid; j++) {
                            int z = ordered_graph[v][j].first;
                            if (set_V.find(z) != set_V.end()) {
                                set_V.erase(z);
                            }
                        }
                    }
                }
                if ((int)S.size() <= K) {
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
        aproximation_server[0] = rand_i() % V;
        int k = 1;
        while(k < K) {
            int next_server = -1, distance = 0;
            for(int v = 0; v < V; v++) {
                int mn = INF;
                bool is_in_the_set = false;
                for(int x: aproximation_server){
                    if (x == -1) break;
                    if (v == x) {
                        is_in_the_set = true;
                        break;
                    }
                    mn = std::min(mn, graph[v][x].second);
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
        if (t >= T) return 0;
        auto &d = dp[t];
        if (d != -1) return d;
        d = 0;
        for (int l = 0; l < L; l++) {
            int lim = std::min(T, t + center_types[l]);
            for (int v = 0; v < V; v++) {
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


    void build() {
        for (int t = 0; t < T; t++) {
            aproximation_server = std::vector<int>(K, -1);
            // dorit(t);
            gon(t);
        }
        for (int k = 0; k < K; k++) {
            for (int &t: dp) t = -1;
            jorge_algos(0);
            int at = 0;
            while(at < T) {
                auto [l, v] = rec[at];
                dna[at].insert({v, l});
                int nxt = std::min(T, at + center_types[l]);
                facility_seg[v].update(at, nxt-1, -1);
                for (int t = at; t < nxt; t++) {
                    vis[t][v] = true;
                }
                at = std::min(T, at + center_types[l]);
            }
        }
        get_fitness();
    }

    int get_fitness() {
        fitness = 0;
        for (int t = 0; t < T; t++) {
            for (auto [v, l]: dna[t]) {
                if (v != -1) {
                    for (int i = t; i < std::min(T, t + center_types[l]); i++) {
                        current_facilities[i][facility_idx[i]++] = v;
                    }
                }
            }
        }

        for (int t = 0; t < T; t++) {
            for (int client: clients[t]) {
                int mn = INF;
                for (int facility: current_facilities[t]) {
                    if (facility == -1) break;
                    mn = std::min(mn, graph[client][facility].second);
                }
                fitness = std::max(mn, fitness);
            }
        }
        
        return fitness;
    }
};


struct Individuo {
    int fitness=-1;
    std::vector<std::set<std::pair<int, int>>> dna, open_facilities;
    std::vector<std::set<int>> used_nodes;
    std::vector<int> mutation_sol, clients_best;
    std::vector<std::pair<int, int>> index_best;

    Individuo() {
        dna = std::vector<std::set<std::pair<int, int>>>(T);
        open_facilities = std::vector<std::set<std::pair<int, int>>>(T);
        used_nodes = std::vector<std::set<int>>(T);
        mutation_sol = std::vector<int>(T, 0);
        clients_best = std::vector<int>(V, INF);
        index_best = std::vector<std::pair<int, int>>(T*K);
    }

    Individuo(std::vector<std::set<std::pair<int, int>>> _matrix, int _fitness) {
        dna = _matrix;
        fitness = _fitness;
        open_facilities = std::vector<std::set<std::pair<int, int>>>(T);
        used_nodes = std::vector<std::set<int>>(T);
        mutation_sol = std::vector<int>(T, 0);
        clients_best = std::vector<int>(V, INF);
        index_best = std::vector<std::pair<int, int>>(T*K);
        count_facilities();
    }

    void count_facilities() {
        for (int t = 0; t < T; t++) open_facilities[t].clear();
        for (int t = 0; t < T; t++) {
            for (auto [v, l]: dna[t]) {
                if (v != -1) {
                    for (int i = t; i < std::min(T, t + center_types[l]); i++) {
                        open_facilities[i].insert({v, t});
                    }
                }
            }
        }
    }

    void mutate() {
        // for (int t = 0; t < T; t++) {
        //     dna[t].clear();
        //     used_nodes[t].clear();
        // }
        // bool ok = true;
        // int t = 0;
        // while (ok) {
        //     int l = rand_i() % L;
        //     ok = false;
        //     for (int i = t; i < T; i++) {
        //         int st = i;
        //         int nd = std::min(T, st + center_types[l]);
        //         int mx = 0;
        //         for (int j = st; j < nd; j++) {
        //             mx = std::max(mx, (int)used_nodes[st].size());
        //         }
        //         if (mx < K) {
        //             std::set<int> available;
        //             for (int v = 0; v < V; v++) {
        //                 bool has = false;
        //                 for (int j = st; j < nd && !has; j++){
        //                     if (used_nodes[j].find(v) != used_nodes[j].end()) {
        //                         has = true;
        //                     }
        //                 }
        //                 if (!has) {
        //                     available.insert(v);
        //                 }
        //             }
        //             if (available.empty()){
        //                 continue;
        //             }
        //             int op = rand_i() % (int)available.size(), v;
        //             for (int j: available) {
        //                 if (!op) {
        //                     v = j;
        //                 }
        //                 op--;
        //             }
        //             for (int u = st; u < nd; u++){
        //                 used_nodes[u].insert(v);
        //             }
        //             dna[i].insert({v, l});
        //             ok = true;
        //             t = i+1;
        //             break;
        //         }
        //         i = nd;
        //     }
        // }

        int k = rand_i() % K + 1;
        for (int t = 0; t < T; t++) {
            open_facilities[t].clear();
            dna[t].clear();
        }
        for (int t = 0; t < T; t++) {
            if ((int)open_facilities[t].size() == K) continue;
            k = rand_i() % (K - (int)open_facilities[t].size()) + 1;
            for (int i = 0; i < k; i++) {
                int l = rand_i() % L;
                int v = rand_i() % V;
                auto it = open_facilities[t].lower_bound({v, -1});
                while (it != open_facilities[t].end() && it->first == v) {
                    v = rand_i() % V;
                    it = open_facilities[t].lower_bound({v, -1});
                }
                int nd = std::min(T, t + center_types[l]);
                for (int st = t; st < nd; st++) {
                    open_facilities[st].insert({v, t});
                }
                dna[t].insert({v, l});
            }
        }

        get_fitness();
        if (give_penalties() != 0) {
            mtx.lock();
            debug(dna);
            debug(open_facilities);
            assert(false);
            mtx.unlock();
        }
    }

    void try_insert(std::set<std::pair<int, int>> gene, Individuo &child, int t) {
        child.count_facilities();
        for (auto [v, l]: gene) {
            int nd = std::min(T, center_types[l] + t);
            bool has = false;
            for (int st = t; st < nd && !has; st++) {
                if (child.open_facilities[st].find({v, -1}) != child.open_facilities[st].end()
                    || (int)child.open_facilities[st].size() == K) {
                    has = true;
                }
            }
            if (has) continue;
            child.dna[t].insert({v, l});
            for (int st = t; st < nd; st++) {
                child.open_facilities[t].insert({v, t});
            }
        }
    }

    void try_insert(std::pair<int, int> alele, Individuo &child, int t) {
        auto [v, l] = alele;
        int nd = std::min(T, center_types[l] + t);
        bool has = false;
        for (int st = t; st < nd && !has; st++) {
            if (child.open_facilities[st].find({v, -1}) != child.open_facilities[st].end()
                || (int)child.open_facilities[st].size() == K) {
                has = true;
            }
        }
        if (has) return;
        child.dna[t].insert({v, l});
        for (int st = t; st < nd; st++) {
            child.open_facilities[t].insert({v, t});
        }
    }

    void cross(Individuo const &p1, Individuo const &p2) {
        Individuo children_1, children_2;
        if (rand() < 0.5) {
            int lim = rand_i() % (T + K);
            for (int t = 0, i = 0; t < T; t++) {
                for (auto alele: p1.dna[t]) {
                    if (i < lim) {
                        try_insert(alele, children_1, t);
                    }
                    else {
                        try_insert(alele, children_2, t);
                    }
                    i++;
                }
                for (auto alele: p2.dna[t]) {
                    if (i < lim) {
                        try_insert(alele, children_2, t);
                    }
                    else {
                        try_insert(alele, children_1, t);
                    }
                    i++;
                }
            }
        }
        else {
            int l = rand_i() % T, r = rand_i() % T;
            for (int t = 0; t < T; t++) {
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
        // bool improve = true;
        // while (improve) {
        //     improve = false;
        //     Individuo save(dna, fitness);
        //     for (int t = 0; t < T; t++) {
        //         int op = rand_i() % 3;
        //         if (!op) {
        //             int v = rand_i() % V;
        //             int l = rand_i() % L;
        //             auto it = dna[t].lower_bound({v, -1});
        //             if (it != dna[t].end() && it->first == v) continue;
        //             dna[t].insert({v, l});
        //         }
        //         else if (!dna[t].empty()) {
        //             int rl, rv;
        //             int k = rand_i() % (int)dna[t].size();
        //             for (auto [i, j]: dna[t]) {
        //                 if (!k) {
        //                     rv = i; rl = j;
        //                     break;
        //                 }
        //                 k--;
        //             }
        //             dna[t].erase(dna[t].lower_bound({rv, rl}));
        //             int v = rand_i() % V;
        //             int l = rand_i() % L;
        //             if (op == 2) {
        //                 auto it = dna[t].lower_bound({v, -1});
        //                 if (it != dna[t].end() && it->first == v) continue;
        //                 dna[t].insert({v, l});
        //             }
        //         }
        //     }
        //     if (get_fitness() < save.fitness) {
        //         improve = true;
        //     }
        //     else {
        //         dna = save.dna;
        //     }
        // }

        Individuo new_sol;
        int cnt = 0;
        for (int t = 0; t < T; t++) {
            cnt += (int)dna[t].size();
        }
        for (int i = 0; i < cnt; i++) index_best[i] = {0, i};

        for (int t = 0; t < T; t++) {
            for (int v = 0; v < V; v++) clients_best[v] = INF;
            for (int client: clients[t]) {
                int best = 0;
                for (int i = 0, c = 0; i < T; i++) {
                    if (dna[i].empty()) {
                        c++; continue;
                    }
                    for (auto [u, l]: dna[i]) {
                        int nd = std::min(T, i + center_types[l]);
                        if (i <= t && t < nd) {
                            if (graph[client][u].second < clients_best[client]) {
                                clients_best[client] = graph[client][u].second;
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
            for (int t = 0, c = 0; t < T; t++) {
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

        for (int t = 0; t < T; t++) {
            used_nodes[t].clear();
        }
        bool ok = true;
        int t = 0;
        while (ok) {
            int l = rand_i() % L;
            ok = false;
            for (int i = t; i < T; i++) {
                int st = i;
                int nd = std::min(T, st + center_types[l]);
                int mx = 0;
                for (int j = st; j < nd; j++) {
                    mx = std::max(mx, (int)used_nodes[st].size());
                }
                if (mx < K) {
                    std::set<int> available;
                    for (int v = 0; v < V; v++) {
                        bool has = false;
                        for (int j = st; j < nd && !has; j++){
                            if (used_nodes[j].find(v) != used_nodes[j].end()) {
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
                        used_nodes[u].insert(v);
                    }
                    new_sol.dna[i].insert({v, l});
                    ok = true;
                    t = i+1;
                    break;
                }
                i = nd;
            }
        }

        for (t = 0; t < T; t++) {
            dna[t] = new_sol.dna[t];
        }

        get_fitness();
    }

    int give_penalties() {
        int total_penalty = 0;

        for (int t = 0; t < T; t++) {
            for (auto [v, l]: dna[t]) {
                if (v != -1) {
                    int nd = std::min(T, t + center_types[l]);
                    for (int st = t; st < nd; st++) {
                        // Check for multifacilities
                        for (auto [u, _]: dna[st]) {
                            if (u == v) {
                                total_penalty += PENALTY;
                            }
                        }
                    }
                    total_penalty -= PENALTY;
                }
            }
        }

        // Check for more than K or no facility
        for (int t = 0; t < T; t++) {
            if (!(int)open_facilities[t].size() || (int)open_facilities[t].size() > K) {
                total_penalty += PENALTY;
            }
        }

        return total_penalty;
    }

    int get_fitness() {
        fitness = 0;
        count_facilities();

        for (int t = 0; t < T; t++) {
            for (int client: clients[t]) {
                int mn = PENALTY;
                for (auto [_, facility]: open_facilities[t]) {
                    if (facility == -1) break;
                    mn = std::min(mn, graph[client][facility].second);
                }
                fitness = std::max(mn, fitness);
            }
        }
        
        fitness += give_penalties();
        
        return fitness;
    }
    
    friend bool operator < (const Individuo lhs, const Individuo rhs){
        return lhs.fitness < rhs.fitness;
    }

    friend std::ostream &operator <<(std::ostream &os, const Individuo &ind) {
        os << "Fitness: " << ind.fitness;
        return os;
    }
};


thread_local std::discrete_distribution<int> d;

/**
 * Main GA struct responsible of running the core of the algorithm.
*/
struct Evolution {
    std::vector<Individuo> population;
    int num_threads, pop_size;
    long long fmax, favg;

    /* Default constructor */
    Evolution(int _pop_size = 100, int _num_threads = 4) {
        num_threads = _num_threads;
        pop_size = _pop_size;
    }

    /**
     * Function to generate a solution inside a chromossome/individuo.
     * @param order Problem specific initialization. (This shall change)
     * @return A single solution.
    */
    InitialSolution generate_individuo() {
        InitialSolution new_ind = InitialSolution();
        new_ind.build();
        return new_ind;
    }

    /**
     * Function to populate a batch of the population.
     * @param l index to start populating.
     * @param r index to end populating.
     * @param order Problem specific initialization. (This shall change)
     * TODO: Use referenced vector instead of struct vector.
    */
    void populate_partially(int l, int r, std::vector<InitialSolution> &solutions) {
        r = std::min(r, pop_size-1);
        for (int i = l; i <= r; i++) {
            solutions[i] = generate_individuo();
        }
    }

    /**
     * Function to initialize the process of population.
    */
    void populate() {
        std::vector<std::thread> threads(num_threads);
        std::vector<InitialSolution> solutions(pop_size);
        int l = 0, step = pop_size/num_threads;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads-1) step = pop_size;
            int r = l + step-1;
            threads[i] = std::thread(&Evolution::populate_partially, this, l, r, std::ref(solutions));
            l += step;
        }
        for (std::thread &th: threads) th.join();
        for (int i = 0; i < pop_size; i++) {
            population[i] = Individuo(solutions[i].dna, solutions[i].fitness);
            debug(solutions[i].fitness);
        }
    }

    /**
     * Select parents on a tournament of 5 candidates.
     * @return A pair o parents.
    */
    static std::pair<Individuo, Individuo> select_parents(std::vector<Individuo> const &this_gen) {
        std::vector<int> parents(5);
        for (int i = 0; i < 5; i++) {
            parents[i] = d(g);
        }
        std::sort(parents.begin(), parents.end(), [&](const int &a, const int &b){
            return this_gen[a] < this_gen[b];
        });
        return {this_gen[parents[0]], this_gen[parents[1]]};
    }

    /**
     * Calculate crossover probability bassed on adaptative calculations.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param fat Fitness of a solution.
     * @return Probability of a crossover to happen.
    */
    static double crossover_prob(int fmax, long long favg, int fat) {
        double k3 = 0.15;
        // return k3;
        double k1 = rand();
        if (k3 > k1) std::swap(k3, k1);
        if (fat < favg) {
            return k3;
        }
        return k1 * std::max((favg - fmax), 1LL) / std::max(fat - fmax, 1);
    }

    /**
     * Calculate mutation probability bassed on adaptative calculations.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param fat Fitness of a solution.
     * @return Probability of a mutation to happen.
    */
    static double mutation_prob(int fmax, long long favg, int fat) {
        double k4 = 0.15;
        // return k4;
        double k2 = rand();
        if (k4 > k2) std::swap(k4, k2);
        if (fat < favg) {
            return k4;
        }
        return k2 * std::max((favg - fmax), 1LL) / std::max(fat - fmax, 1);
    }

    /**
     * Function to crossover a batch of the population.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param l index to start the crossover.
     * @param r index to end the crossover.
     * @param next_gen Next generation of solutions.
    */
    static void run_crossover_partially(int fmax, long long favg, int l, int r, std::vector<Individuo> &next_gen, std::vector<Individuo> const &this_gen) {
        for (int i = l; i <= r; i++) {
            if (rand() < crossover_prob(fmax, favg, this_gen[i].fitness)) {
                std::pair<Individuo, Individuo> parents = select_parents(this_gen);
                next_gen[i].cross(parents.first, parents.second);
            }
        }
    }

    void run_crossover(int elitism_size, std::vector<Individuo> &next_gen) {
        int l = elitism_size, step = (pop_size - elitism_size) / num_threads;
        std::vector<std::thread> threads(num_threads);
        for (int i = 0; i < num_threads; i++) {
            int r;
            if (i == num_threads - 1) r = pop_size-1;
            else r = l + step - 1;
            threads[i] = std::thread(run_crossover_partially, fmax, favg, l, r, std::ref(next_gen), std::cref(population));
            l += step;
        }
        for (std::thread &th: threads) th.join();
    }

    /**
     * Function to mutation a batch of the population.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param l index to start the mutation.
     * @param r index to end the mutation.
     * @param next_gen Next generation of solutions.
    */
    static void run_mutation_partially(int fmax, long long favg, int l, int r, std::vector<Individuo> &next_gen, std::vector<Individuo> const &this_gen) {
        for (int i = l; i <= r; i++) {
            if (rand() < mutation_prob(fmax, favg, this_gen[i].fitness)) {
                next_gen[i].mutate();
            }
        }
    }


    void run_mutation(int elitism_size, std::vector<Individuo> &next_gen) {
        int l = elitism_size, step = (pop_size - elitism_size) / num_threads;
        std::vector<std::thread> threads(num_threads);
        for (int i = 0; i < num_threads; i++) {
            int r;
            if (i == num_threads - 1) r = pop_size-1;
            else r = l + step - 1;
            threads[i] = std::thread(run_mutation_partially, fmax, favg, l, r, std::ref(next_gen), std::cref(population));
            l += step;
        }
        for (std::thread &th: threads) th.join();
    }


    /**
     * Function to run a local search on a batch of the population.
     * @param l index to start the local search.
     * @param r index to end the local search.
     * @param next_gen Next generation of solutions.
    */
    static void run_local_search_partially(int l, int r, std::vector<Individuo> &next_gen) {
        for (int i = l; i <= r; i++) {
            next_gen[i].local_search();
        }
    }


    void run_local_search(int elitism_size, std::vector<Individuo> &next_gen) {
        int l = elitism_size, step = (pop_size - elitism_size) / num_threads;
        std::vector<std::thread> threads(num_threads);
        for (int i = 0; i < num_threads; i++) {
            int r;
            if (i == num_threads - 1) r = pop_size-1;
            else r = l + step - 1;
            threads[i] = std::thread(run_local_search_partially, l, r, std::ref(next_gen));
            l += step;
        }
        for (std::thread &th: threads) th.join();
    }


    /**
     * Function to handle the evolutionary process.
    */
    std::pair<Individuo, double> run_evo(
        int fitness_limit = -1, double elitism = 0.1, double time_limit = 3600, bool verbose = false
    ) {
        population.assign(pop_size, Individuo());
        int elitism_size = std::min(int(pop_size * elitism), 3);
        double start = clock();
        debug("creatining population");
        populate();

        int best = INF;
        int gen = 0;
        std::set<std::vector<std::set<std::pair<int, int>>>> past_solutions;
        for (Individuo &ind: population) past_solutions.insert(ind.dna);
        double nd = (clock() - start) / CLOCKS_PER_SEC;
        while ((double)(clock() - start) / CLOCKS_PER_SEC < time_limit) {
            std::vector<Individuo> next_gen(pop_size);
            debug(gen);
            std::sort(population.begin(), population.end());
            best = std::min(best, population[0].fitness);
            if (best < population[0].fitness) {
                nd = (clock() - start) / CLOCKS_PER_SEC;
                best = population[0].fitness;
            }

            if (best == fitness_limit) {
                break;
            }

            assert((int)population.size() == pop_size);
            int worst = population.back().fitness;

            fmax = INF;
            favg = 0;
            int sample = std::min(pop_size, 10);
            for (int i = 0; i < sample; i++) {
                Individuo ind = population[i];
                fmax = std::min(fmax, 1ll * ind.fitness);
                favg += 1ll * ind.fitness;
                debug(ind);
            }
            favg /= sample;

            if (verbose) {
                std::map<std::vector<std::set<std::pair<int, int>>>, int> occ;
                int ns = 0;
                for (Individuo &ind: population) {
                    if (past_solutions.find(ind.dna) == past_solutions.end()) ns++;
                    occ[ind.dna]++;
                }

                double gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2.0) * 100;
                std::cout << std::string(100, ' ') << '\r';
                std::cout.flush();
                std::cout << "Generation " << gen << ", Fitness: min=" << fmax << ", mean=" << favg
                    << ", max=" << worst << ", Gap:" << std::fixed << std::setprecision(6) << gap  << ", dif sol: " << (int)occ.size() << ", new sol: " << ns << "\r";
                std::cout.flush();

                past_solutions.clear();
                for (Individuo &ind: population) {
                    past_solutions.insert(ind.dna);
                }
            }
            std::vector<int> weigths(pop_size);
            int cur_weight = pop_size + 5;
            int ant = best;
            for (int i = 0; i < pop_size; i++) {
                if (population[i].fitness > ant) {
                    ant = population[i].fitness;
                    cur_weight--;
                }
                weigths[i] = cur_weight;
            }
            d = std::discrete_distribution<int>(weigths.begin(), weigths.end());

            for (int i = 0; i < pop_size; i++) {
                next_gen[i] = population[i];
            }
            debug("starting local search");
            run_local_search(0, next_gen);

            debug("starting crossover");
            run_crossover(elitism_size, next_gen);

            
            debug("starting mutation");
            run_mutation(elitism_size, next_gen);


            debug("moving");
            std::sort(next_gen.begin(), next_gen.end());
            for (int i = elitism_size; i < pop_size; i++) {
                population[i] = next_gen[i-elitism_size];
            }
            gen++;
        }

        if (verbose) {
            std::cout << std::string(90, ' ') << '\r';
            std::cout.flush();
            double gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2.0) * 100;
            std::cout << "Generation " << gen << ", Fitness =" << best << std::fixed << std::setprecision(6) << ", Gap: " << gap << "\n";
        }

        return {population[0], nd};
    }
};

// 1: path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 3);

    std::string dir = argv[1];
    read(dir + argv[2]);


    int fitness_limit = 125;
    double elitism = .1;
    int pop_size = 1000;
    int num_threads = 10;
    double time_limit = 60*5;
    bool verbose = true;
    Evolution ga(pop_size, num_threads);

    auto [ind, t] = ga.run_evo(fitness_limit, elitism, time_limit, verbose);
    std::cout << ind << '\n';
    std::cout << std::fixed << std::setprecision(5) << t << '\n';




    return 0;
}
