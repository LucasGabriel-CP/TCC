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
#include <unordered_set>
#include <unordered_map>
#include "util/Constants.cpp"
#include "util/DebugTemplate.cpp"

int V, T, L, K, E;
std::vector<std::vector<std::pair<int, int>>> graph;      // dist(i, j)
std::vector<std::vector<int>> clients;          // Clientes no instante t
std::vector<int> center_types;                  // Tipos de centro

void read(std::string const &instance_path) {
    std::ifstream in;

    debug(instance_path);
    in.open(instance_path);

    in >> V >> T >> L >> K;

    graph.assign(V, std::vector<std::pair<int, int>>(V));
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            graph[i][j].first = j;
            in >> graph[i][j].second;
        }
        debug(graph[i]);
        // std::sort(graph[i].begin(), graph[i].end(), [](auto a, auto b){
        //     return a.second < b.second;
        // });
    }

    center_types.resize(L);
    for (int i = 0; i < L; i++) {
        in >> center_types[i];
    }

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
    std::vector<std::vector<int>> matrix, current_facilities; // L por T
    std::vector<std::pair<int, int>> rec;
    std::vector<std::vector<bool>> vis;
    std::vector<int> gon_server, dp;
    std::vector<SegTree> facility_seg;


    InitialSolution() {
        matrix = std::vector<std::vector<int>>(L, std::vector<int>(T, -1));
        current_facilities = std::vector<std::vector<int>>(T, std::vector<int>(L, -1));
        rec = std::vector<std::pair<int, int>>(T);
        gon_server = std::vector<int>(K, -1);
        dp = std::vector<int>(T, -1);
        vis = std::vector<std::vector<bool>>(T, std::vector<bool>(V, false));
        facility_seg = std::vector<SegTree>(V, SegTree(T));
    }

    void gon(int t) {
        gon_server[0] = rand_i() % V;
        int k = 1;
        debug(graph);
        while(k < K) {
            int next_server = -1, distance = 0;
            debug(gon_server);
            for(int v = 0; v < V; v++) {
                int mn = INF;
                bool is_in_the_set = false;
                for(int x: gon_server){
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
            gon_server[k++] = next_server;
        }
        for (int i: gon_server) {
            debug(i);
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
            gon_server = std::vector<int>(K, -1);
            gon(t);
        }
        for (int k = 0; k < K; k++) {
            for (int &t: dp) t = -1;
            jorge_algos(0);
            int at = 0;
            while(at < T) {
                auto [l, v] = rec[at];
                matrix[l][at] = v;
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
            for (int l = 0; l < L; l++) {
                if (matrix[l][t] != -1) {
                    for (int i = t; i < std::min(T, t + center_types[l]); i++) {
                        current_facilities[i][l] = matrix[l][t];
                    }
                }
            }
        }

        for (int t = 0; t < T; t++) {
            for (int client: clients[t]) {
                int mn = INF;
                for (int facility: current_facilities[t]) {
                    if (facility == -1) continue;
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
    std::vector<std::vector<int>> matrix, current_facilities;
    std::vector<std::vector<bool>> open_facilities, used_nodes; // L por T
    std::vector<int> open_time;


    Individuo() {
        matrix = std::vector<std::vector<int>>(L, std::vector<int>(T, -1));
        current_facilities = std::vector<std::vector<int>>(T, std::vector<int>(L, -1));
        open_facilities = std::vector<std::vector<bool>>(L, std::vector<bool>(T, false));
        used_nodes = std::vector<std::vector<bool>>(V, std::vector<bool>(T, false));
        open_time = std::vector<int>(T, 0);
    }

    Individuo(std::vector<std::vector<int>> _matrix, int _fitness) {
        matrix = _matrix;
        fitness = _fitness;
        open_time = std::vector<int>(T, 0);
    }

    void fix_this() {
        for (int l = 0; l < L; l++) {
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    int lim = std::min(T, t + center_types[t]);
                    for (int i = t+1; i < lim; i++) {
                        matrix[l][t] = -1;
                    }
                }
            }
        }
    }

    void mutate() {
        if (rand() < 0.5) {
            int idx = rand_i() % L;
            int l = rand_i() % T, r = rand_i() % T;
            if (l > r) std::swap(l, r);

            for (; l < r; l++, r--) {
                std::swap(matrix[idx][l], matrix[idx][r]);
            }
        }
        else {
            int idx_1 = rand_i() % L, idx_2 = rand_i() % L;
            int l = rand_i() % T, r = rand_i() % T;
            if (l > r) std::swap(l, r);
            for (; l <= r; l++) {
                std::swap(matrix[idx_1][l], matrix[idx_2][l]);
            }
        }

        get_fitness();
    }

    static Individuo cross(Individuo const &p1, Individuo const &p2) {
        Individuo children_1, children_2;
        if (rand() < 0.5) {
            int lim = 0;
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < T; j++) {
                    if (j > lim) {
                        children_1.matrix[i][j] = p1.matrix[i][j];
                        children_2.matrix[i][j] = p2.matrix[i][j];
                    }
                    else {
                        children_1.matrix[i][j] = p2.matrix[i][j];
                        children_2.matrix[i][j] = p1.matrix[i][j];
                    }
                }
                lim++;
            }
        }
        else {
            for (int i = 0; i < L; i++) {
                if (rand() < 0.5) {
                    children_1.matrix[i] = p1.matrix[i];
                    children_2.matrix[i] = p2.matrix[i];
                }
                else {
                    children_1.matrix[i] = p2.matrix[i];
                    children_2.matrix[i] = p1.matrix[i];
                }
            }
        }

        // Select a random child to return
        if (rand() < 0.5) {
            children_1.fix_this();
            children_1.get_fitness();
            return children_1;
        }
        children_2.fix_this();
        children_2.get_fitness();
        return children_2;
    }

    bool validate() {        
        #ifdef PIZZA
            // Verify solution
        #endif

        return true;
    }

    int give_penalties() {
        int total_penalty = 0;

        return total_penalty;
    }

    void local_search() {
        for (int &i: open_time) i = 0;
        for (int l = 0; l < L; l++){
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    int lim = std::min(T, t + center_types[l]);
                    for (int i = t; i < lim; i++) {
                        open_facilities[l][i] = true;
                        used_nodes[matrix[l][t]][i] = true;
                    }
                    open_time[t]++;
                }
            }
        }

        for (int t = 0; t < T; t++) {
            if (open_time[t] < K) {
                int l = rand_i() % L;
                int v = rand_i() % V;
                int lhs = std::min(T, t + center_types[l]) - 1;
                int rhs = std::min(0, t - center_types[l] + 1);
                int i;
                for (i = lhs; i <= rhs && !open_facilities[l][i] && !used_nodes[v][i]; i++);
                if (i > rhs) {
                    matrix[l][t] = v;
                    t = rhs;
                }
            }
        }
    }

    int get_fitness() {
        fitness = 0;
        current_facilities = std::vector<std::vector<int>>(T, std::vector<int>(L, -1));
        for (int t = 0; t < T; t++) {
            int idx = 0;
            for (int l = 0; l < L; l++) {
                if (matrix[l][t] != -1) {
                    for (int i = t; i < std::min(T, t + center_types[l]); i++) {
                        current_facilities[t][idx++] = matrix[l][t];
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
            if (i == num_threads-1) step = pop_size; // Give remainder to last thread, prob better to give this to the first one.
            int r = l + step-1;
            threads[i] = std::thread(&Evolution::populate_partially, this, l, r, std::ref(solutions));
            l += step;
        }
        for (std::thread &th: threads) th.join();
        // for (int i = 0; i < pop_size; i++) {
        //     solutions[i] = generate_individuo();
        // }
        for (int i = 0; i < pop_size; i++) {
            population[i] = Individuo(solutions[i].matrix, solutions[i].fitness);
            debug(solutions[i].fitness);
        }
    }

    /**
     * Select parents on a tournament of 5 candidates.
     * @return A pair o parents.
    */
    std::pair<Individuo, Individuo> select_parents() {
        std::vector<int> parents(5);
        for (int i = 0; i < 5; i++) {
            parents[i] = d(g);
        }
        std::sort(parents.begin(), parents.end(), [&](const int &a, const int &b){
            return population[a] < population[b];
        });
        return {population[parents[0]], population[parents[1]]};
    }

    /**
     * Calculate crossover probability bassed on adaptative calculations.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param fat Fitness of a solution.
     * @return Probability of a crossover to happen.
    */
    double crossover_prob(double fmax, double favg, double fat) {
        double k3 = 0.15;
        double k1 = rand();
        if (fat < favg) {
            return k3;
        }
        return k1 * std::max((favg - fmax), 1.0) / std::max(fat - fmax, 1.0);
    }

    /**
     * Calculate mutation probability bassed on adaptative calculations.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param fat Fitness of a solution.
     * @return Probability of a mutation to happen.
    */
    double mutation_prob(double fmax, double favg, double fat) {
        double k4 = 0.15;
        double k2 = rand();
        if (fat < favg) {
            return k4;
        }
        return k2 * std::max((favg - fmax), 1.0) / std::max(fat - fmax, 1.0);
    }

    /**
     * Function to crossover a batch of the population.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param l index to start the crossover.
     * @param r index to end the crossover.
     * @param next_gen Next generation of solutions.
    */
    void run_crossover_partially(double fmax, double favg, int l, int r, std::vector<Individuo> &next_gen) {
        for (int i = l; i <= r; i++) {
            if (rand() < crossover_prob(fmax, favg, population[i].fitness)) {
                std::pair<Individuo, Individuo> parents = select_parents();
                Individuo child = Individuo::cross(parents.first, parents.second);
                next_gen[i] = child;
            }
            else {
                next_gen[i] = population[i];
            }
        }
    }

    /**
     * Function to mutation a batch of the population.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param l index to start the mutation.
     * @param r index to end the mutation.
     * @param next_gen Next generation of solutions.
    */
    void run_mutation_partially(double fmax, double favg, int l, int r, std::vector<Individuo> &next_gen) {
        for (int i = l; i <= r; i++) {
            if (rand() < mutation_prob(fmax, favg, population[i].fitness)) {
                next_gen[i].mutate();
            }
        }
    }


    /**
     * Function to run a local search on a batch of the population.
     * @param l index to start the local search.
     * @param r index to end the local search.
     * @param next_gen Next generation of solutions.
    */
    void run_local_search_partially(int l, int r, std::vector<Individuo> &next_gen) {
        for (int i = l; i <= r; i++) {
            if (rand() < 0.5) {
                next_gen[i].local_search();
            }
        }
    }


    /**
     * Function to handle the evolutionary process.
    */
    std::pair<Individuo, int> run_evo(
        int fitness_limit = -1, double elitism = 0.1, double time_limit = 3600, bool verbose = false
    ) {
        population.assign(pop_size, Individuo());
        int elitism_size = std::min(int(pop_size * elitism), 3);
        int start = clock();
        debug("creatining population");
        populate();
        debug("population created");

        int best = INF;
        int gen = 0;
        while ((double)(clock() - start) / CLOCKS_PER_SEC < time_limit) {
            std::sort(population.begin(), population.end());
            best = std::min(best, population[0].fitness);

            if (best == fitness_limit) {
                break;
            }

            assert((int)population.size() == pop_size);
            double worst = population.back().fitness;

            int fmax = INF, favg = 0;
            for (Individuo &ind: population) {
                fmax = std::min(fmax, ind.fitness);
                favg += ind.fitness;
            }
            favg /= pop_size;

            if (verbose) {
                double gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2.0) * 100;
                std::cout << std::string(80, ' ') << '\r';
                std::cout.flush();
                std::cout << "Generation " << gen << ", Fitness: min=" << fmax << ", mean=" << favg
                    << ", max=" << worst << ", Gap:" << std::fixed << std::setprecision(6) << gap << "\r";
                std::cout.flush();
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

            std::vector<std::thread> threads(num_threads);
            std::vector<Individuo> next_gen(pop_size);
            int l = elitism_size, step = (pop_size - elitism_size) / num_threads;
            int other_step = (pop_size - elitism_size) % num_threads;
            debug("starting crossover");
            for (int i = 0; i < num_threads; i++) {
                int r;
                if (!i) r = l + other_step + step - 1;
                else r = l + step - 1;
                threads[i] = std::thread(&Evolution::run_crossover_partially, this, fmax, favg, l, r, std::ref(next_gen));
                l += step;
                if (!i) l += other_step;
            }
            for (std::thread &th: threads) th.join();
            
            debug("starting mutation");
            l = elitism_size;
            for (int i = 0; i < num_threads; i++) {
                int r;
                if (!i) r = l + other_step + step - 1;
                else r = l + step - 1;
                threads[i] = std::thread(&Evolution::run_mutation_partially, this, fmax, favg, l, r, std::ref(next_gen));
                l += step;
                if (!i) l += other_step;
            }
            for (std::thread &th: threads) th.join();

            l = elitism_size;
            for (int i = 0; i < num_threads; i++) {
                int r;
                if (!i) r = l + other_step + step - 1;
                else r = l + step - 1;
                threads[i] = std::thread(&Evolution::run_local_search_partially, this, l, r, std::ref(next_gen));
                l += step;
                if (!i) l += other_step;
            }
            for (std::thread &th: threads) th.join();

            for (int i = 0; i < elitism_size; i++) next_gen[i] = population[i];
            population = std::move(next_gen);

            gen++;
        }

        debug("batat");
        if (verbose) {
            double gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2.0) * 100;
            std::cout << "Generation " << gen << ", Fitness =" << best << std::fixed << std::setprecision(6) << ", Gap: " << gap << "\n";
        }

        int nd = (clock() - start) / CLOCKS_PER_SEC;

        return {population[0], nd};
    }
};


// 1: path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 3);

    std::string dir = argv[1];
    read(dir + argv[2]);


    int fitness_limit = 20;
    double elitism = .2;
    int pop_size = 24;
    int num_threads = 6;
    double time_limit = 3600;
    bool verbose = true;
    Evolution ga(pop_size, num_threads);

    auto [ind, val] = ga.run_evo(fitness_limit, elitism, time_limit, verbose);
    std::cout << ind << '\n';




    return 0;
}
