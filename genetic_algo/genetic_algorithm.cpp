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

thread_local int V, T, L, K, E;
thread_local std::vector<std::vector<std::pair<int, int>>> graph;      // dist(i, j)
thread_local std::vector<std::vector<int>> clients;          // Clientes no instante t
thread_local std::vector<int> center_types;                  // Tipos de centro

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
        std::sort(graph[i].begin(), graph[i].end(), [](auto a, auto b){
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
        debug(clients[t]);
    }

    in.close();
}



struct Segtree {
	std::vector<long long> seg, lazy;
	int n, LOG;

    Segtree(int _n) {
        n = _n;
        seg.assign(n, 0ll);
        lazy.assign(n, 0ll);
        LOG = std::ceil(log2(n));
    }

	long long junta(long long a, long long b) {
		return a+b;
	}

	void push(int p, long long x, int tam, bool prop=1) {
		seg[p] += x*tam;
		if (prop and p < n) lazy[p] += x;
	}

	void up(int p) {
		for (int tam = 2; p /= 2; tam *= 2) {
			seg[p] = junta(seg[2*p], seg[2*p+1]);
			push(p, lazy[p], tam, 0);
		}
	}

	void prop(int p) {
		int tam = 1 << (LOG-1);
		for (int s = LOG; s; s--, tam /= 2) {
			int i = p >> s;
			if (lazy[i]) {
				push(2*i, lazy[i], tam);
				push(2*i+1, lazy[i], tam);
				lazy[i] = 0;
			}
		}
	}

	long long query(int a, int b) {
		long long ret = 0;
		for (prop(a+=n), prop(b+=n); a <= b; ++a/=2, --b/=2) {
			if (a%2 == 1) ret = junta(ret, seg[a]);
			if (b%2 == 0) ret = junta(ret, seg[b]);
		}
		return ret;
	}

	void update(int a, int b, int x) {
		int a2 = a += n, b2 = b += n, tam = 1;
		for (; a <= b; ++a/=2, --b/=2, tam *= 2) {
			if (a%2 == 1) push(a, x, tam);
			if (b%2 == 0) push(b, x, tam);
		}
		up(a2), up(b2);
	}
};


struct Individuo {
    int fitness=-1;
    std::vector<std::vector<int>> matrix; // L por T

    Individuo() {
        matrix.assign(L, std::vector<int>(T, -1));
    }

    void dorit(int t, std::vector<Segtree> &facility_seg) {
        int low = 0, high = V-1;
        std::vector<int> ans;
        while (low < high) {
            int mid = (low + high) / 2;
            std::vector<int> S;
            std::unordered_set<int> set_V;
            for (int i = 0; i < V; i++) set_V.insert(i);
            while (!set_V.empty()) {
                int x = *set_V.begin();
                S.push_back(x);
                for (int i = 0; i <= mid; i++) {
                    int v = graph[x][i].first;
                    if (set_V.find(v) != set_V.end()) {
                        set_V.erase(v);
                        for (int j = 0; j <= mid; j++) {
                            int z = graph[v][j].first;
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


    int jorge_algos(
        int t, std::vector<int> &dp, std::vector<Segtree> &facility_seg,
        std::vector<std::pair<int, int>> &rec, std::vector<std::vector<bool>> &vis
    ) {
        if (t >= T) return 0;
        auto &d = dp[t];
        if (d != -1) return d;
        d = 0;
        for (int l = 0; l < L; l++) {
            int lim = std::min(T, t + center_types[l]);
            for (int v = 0; v < V; v++) {
                if (vis[t][v]) continue;
                int qnt = facility_seg[v].query(t, lim-1);
                int aux = jorge_algos(lim, dp, facility_seg, rec, vis) + qnt;
                if (aux > d)  {
                    rec[t] = {l, v};
                    d = aux;
                }
            }
        }
        return d;
    }


    void build() {
        std::vector<Segtree> facility_seg(V, Segtree(T));
        debug(V, T, K);
        for (int t = 0; t < T; t++) {
            dorit(t, facility_seg);
        }

        std::vector<std::pair<int, int>> rec(T);
        std::vector<std::vector<bool>> vis(T, std::vector<bool>(V, false));
        for (int k = 0; k < K; k++) {
            std::vector<int> dp(T, -1);
            jorge_algos(0, dp, facility_seg, rec, vis);
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
        
        // assert(validate());
        get_fitness();
        mtx.lock();
        debug(fitness);
        mtx.unlock();
    }

    void fix_this() {
        for (int l = 0; l < L; l++) {
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    int lim = std::min(T, t + center_types[t]);
                    for (int i = t; i < lim; i++) {
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
        std::vector<Segtree> open_facilities(L, Segtree(T)), used_nodes(V, Segtree(T));
        std::vector<int> open_time(T, 0);
        for (int l = 0; l < L; l++){
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    int lim = std::min(T, t + center_types[l]);
                    open_facilities[l].update(t, lim-1, 1);
                    used_nodes[matrix[l][t]].update(t, lim-1, 1);
                    open_time[t]++;
                }
            }
        }

        for (int t = 0; t < T; t++) {
            if (open_time[t] < K) {
                int l = rand_i() % L;
                int lhs = std::min(T, t + center_types[l]) - 1;
                int rhs = std::min(0, t - center_types[l] + 1);
                int qnt = open_facilities[l].query(lhs, rhs);
                if (!qnt) {
                    int v = rand_i() % V;
                    qnt = used_nodes[v].query(lhs, rhs);
                    if (!qnt) {
                        matrix[l][t] = v;
                        t = rhs;
                    }
                }
            }
        }
    }

    int get_fitness() {
        fitness = 0;
        std::vector<std::vector<int>> open_facilities(T);
        for (int l = 0; l < L; l++) {
            for (int t = 0; t < T; t++) {
                if (matrix[l][t] != -1) {
                    for (int i = t; i < std::min(T, t + center_types[l]); i++) {
                        open_facilities[t].push_back(matrix[l][t]);
                    }
                }
            }
        }

        for (int t = 0; t < T; t++) {
            for (int client: clients[t]) {
                int mn = INF;
                for (int facility: open_facilities[t]) {
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
    int num_t;

    /* Default constructor */
    Evolution() { }

    /**
     * Problem Specific Constructor
     * @param _t Size of the chromossome
    */
    Evolution(int _t) {
        num_t = _t;
    }

    /**
     * Function to generate a solution inside a chromossome/individuo.
     * @param order Problem specific initialization. (This shall change)
     * @return A single solution.
    */
    Individuo generate_individuo() {
        Individuo new_ind = Individuo();
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
    void populate_partially(int l, int r) {
        r = std::min(r, (int)population.size()-1);
        for (int i = l; i <= r; i++) {
            population[i] = generate_individuo();
        }
    }

    /**
     * Function to initialize the process of population.
     * @param num_threads Number of threads to be created.
    */
    void populate(int num_threads) {
        std::vector<std::thread> threads(num_threads);
        // std::vector std::vector<Individuo> initial_solution((int)population.size());
        int l = 0, step = (int)population.size()/num_threads;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads-1) step = (int)population.size(); // Give remainder to last thread, prob better to give this to the first one.
            int r = l + step-1;
            threads[i] = std::thread(&Evolution::populate_partially, this, l, r);
            l += step;
        }
        for (std::thread &th: threads) th.join();
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
        int fitness_limit = -1, double elitism = 0.1, int pop_size = 100, int num_threads = 4,
        double time_limit = 3600, bool verbose = false
    ) {
        population.resize(pop_size);
        int elitism_size = std::min(int(pop_size * elitism), 3);
        int start = clock();
        debug("creatining population");
        populate(num_threads);
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
            for (int i = 0; i < num_threads; i++) {
                int r;
                if (!i) r = l + other_step + step - 1;
                else r = l + step - 1;
                threads[i] = std::thread(&Evolution::run_crossover_partially, this, fmax, favg, l, r, std::ref(next_gen));
                l += step;
                if (!i) l += other_step;
            }
            for (std::thread &th: threads) th.join();
            
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


        if (verbose) {
            double gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2.0) * 100;
            std::cout << "Generation " << gen << ", Fitness =" << best << std::fixed << std::setprecision(6) << ", Gap: " << gap << "\n";
        }

        int nd = (clock() - start) / CLOCKS_PER_SEC;

        return {population[0], nd};
    }
};
