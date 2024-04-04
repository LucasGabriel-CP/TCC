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
#include "individuo.cpp"
#include "initial_solution.cpp"

thread_local std::discrete_distribution<int> d;

/**
 * Main GA struct responsible of running the core of the algorithm.
*/
struct Evolution {
    std::vector<Individuo> population;
    LeasingProblem problem;
    int num_threads, pop_size;
    long long fmax, favg;

    /* Default constructor */
    Evolution(LeasingProblem const &_problem, int _pop_size = 100, int _num_threads = 4) {
        problem = _problem;
        num_threads = _num_threads;
        pop_size = _pop_size;
    }

    /**
     * Function to generate a solution inside a chromossome/individuo.
     * @param order Problem specific initialization. (This shall change)
     * @return A single solution.
    */
    InitialSolution generate_individuo() {
        InitialSolution new_ind = InitialSolution(problem);
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
        std::vector<InitialSolution> solutions(pop_size, InitialSolution(problem));
        int l = 0, step = pop_size/num_threads;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads-1) step = pop_size;
            int r = l + step-1;
            threads[i] = std::thread(&Evolution::populate_partially, this, l, r, std::ref(solutions));
            l += step;
        }
        for (std::thread &th: threads) th.join();
        for (int i = 0; i < pop_size; i++) {
            population[i] = Individuo(solutions[i].dna, solutions[i].fitness, problem);
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
    static double crossover_prob(long long fmax, long long favg, int fat) {
        double k3 = 0.15;
        // return k3;
        double k1 = rand();
        if (k3 > k1) std::swap(k3, k1);
        if (fat < favg) {
            return k3;
        }
        return k1 * std::max((favg - fmax), 1LL) / std::max(fat - fmax, 1LL);
    }

    /**
     * Calculate mutation probability bassed on adaptative calculations.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param fat Fitness of a solution.
     * @return Probability of a mutation to happen.
    */
    static double mutation_prob(long long fmax, long long favg, int fat) {
        double k4 = 0.15;
        // return k4;
        double k2 = rand();
        if (k4 > k2) std::swap(k4, k2);
        if (fat < favg) {
            return k4;
        }
        return k2 * std::max((favg - fmax), 1LL) / std::max(fat - fmax, 1LL);
    }

    /**
     * Function to crossover a batch of the population.
     * @param fmax Best fitness.
     * @param favg Average fitness.
     * @param l index to start the crossover.
     * @param r index to end the crossover.
     * @param next_gen Next generation of solutions.
    */
    static void run_crossover_partially(long long fmax, long long favg, int l, int r, std::vector<Individuo> &next_gen, std::vector<Individuo> const &this_gen) {
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
    static void run_mutation_partially(long long fmax, long long favg, int l, int r, std::vector<Individuo> &next_gen, std::vector<Individuo> const &this_gen) {
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
            next_gen[i].local_searchv2();
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
        population.assign(pop_size, Individuo(problem));
        int elitism_size = std::min(int(pop_size * elitism), 3);
        debug("creatining population");
        populate();

        long long best = INF;
        int gen = 0;
        time_t start, nd;
        time(&start);
        time(&nd);
        std::vector<Individuo> next_gen(pop_size);
        time_t at;
        time(&at);
        while ((double)(at - start) < time_limit) {
            debug(gen);
            std::sort(population.begin(), population.end());
            if (best > population[0].fitness) {
                time(&nd);
                best = population[0].fitness;
            }

            if (best == fitness_limit) {
                break;
            }

            assert((int)population.size() == pop_size);
            long long worst = population.back().fitness;

            fmax = INF;
            favg = 0;
            for (int i = 0; i < pop_size; i++) {
                Individuo ind = population[i];
                fmax = std::min(fmax, 1ll * ind.fitness);
                favg += 1ll * ind.fitness;
                // debug(ind);
            }
            favg /= pop_size;

            if (verbose) {
                double gap = double(population[0].fitness - fitness_limit) / fitness_limit * 100;
                std::cout << std::string(100, ' ') << '\r';
                std::cout.flush();
                std::cout << "Generation " << gen << ", Fitness: min=" << fmax << ", mean=" << favg
                    << ", max=" << worst << ", Gap:" << std::fixed << std::setprecision(6) << gap  << "\r";
                std::cout.flush();
            }

            std::vector<int> weigths(pop_size);
            int cur_weight = pop_size + 5;
            long long ant = best;
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


            debug("sorting", next_gen.size());
            std::sort(next_gen.begin(), next_gen.end());
            debug("moving");
            for (int i = elitism_size; i < pop_size; i++) {
                population[i] = next_gen[i-elitism_size];
            }
            debug("moved");
            gen++;
            time(&at);
        }

        if (true) {
            std::cout << std::string(90, ' ') << '\r';
            std::cout.flush();
            double gap = double(population[0].fitness - fitness_limit) / fitness_limit * 100;
            std::cout << "Generation " << gen << ", Fitness =" << best << std::fixed << std::setprecision(6) << ", Gap: " << gap << "\n";
        }

        assert(!population[0].give_penalties());

        return {population[0], (double)(nd - start)};
    }
};