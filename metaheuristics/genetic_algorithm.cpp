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
#include "neighbor.cpp"
#include "variable_neighborhood.cpp"

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
        double k3 = 0.25;
        // return k3;
        double k1 = get_rand(0.0, k3);
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
        double k4 = 0.25;
        // return k4;
        double k2 = get_rand(0.0, k4);
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
            if (rand() < 0.5) {
                next_gen[i].first_improve_ls();
            }
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
        long long fitness_limit = -1, double elitism = 0.1, double time_limit = 3600, bool verbose = false
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
        std::vector<Individuo> next_gen(pop_size, Individuo(problem));
        time_t at;
        time(&at);
        int cnt = 0;
        while ((double)(at - start) < time_limit) {
            debug(gen);
            std::sort(population.begin(), population.end());
            if (best > population[0].fitness) {
                time(&nd);
                best = population[0].fitness;
                cnt = 0;
            }
            else {
                cnt++;
            }
            if (cnt == 3000) {
                for (int i = 0; i < elitism_size; i++) {
                    next_gen[i] = population[i];
                }
                populate();
                for (int i = 0; i < elitism_size; i++) {
                    population[i] = next_gen[i];
                }
                cnt = 0;
            }

            double gap = double(population[0].fitness - fitness_limit) / population[0].fitness * 100;
            if (gap < 0.01) {
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

        population[0].valid();

        return {population[0], (double)(nd - start)};
    }

    static void run_vns(
        std::vector<Individuo> const &population, LeasingProblem const &problem, std::vector<Neighbor> &vns_solutions,
        int time_limit, long long fitness_limit, int num_threads, int pop_size
    ) {
        for (int i = 0; i < pop_size; i++) {
            vns_solutions[i].dna = population[i].dna;
            vns_solutions[i].fitness = population[i].fitness;
        }

        auto runner = [&](int _l, int _r) {
            VNS model(problem);
            model.run_multisol(vns_solutions, _l, _r, time_limit, fitness_limit);
        };
        int l = 0, divs = pop_size / num_threads;
        std::vector<std::thread> threads(num_threads);
        for (int i = 0; i < num_threads; i++) {
            int r = std::min(pop_size-1, l + divs - 1);
            if (i < pop_size % num_threads) r++;
            threads[i] = std::thread(runner, l, r);
            l += divs;
            if (i < pop_size % num_threads) l++;
        }
        for (std::thread &th: threads) th.join();
        std::sort(vns_solutions.begin(), vns_solutions.end());
    }

    std::pair<Individuo, double> run_hibrid(
        long long fitness_limit = -1, double elitism = 0.1, double main_tl = 3600, double second_tl = 360, bool verbose = false
    ) {
        auto run_ga = [&](int elitism_size) {
            time_t start, nd, at;
            time(&start);
            time(&nd);
            time(&at);
            std::vector<Individuo> next_gen(pop_size, Individuo(problem));
            long long best = INF;

            while ((double)(at - start) < second_tl) {
                std::sort(population.begin(), population.end());
                if (best > population[0].fitness) {
                    time(&nd);
                    best = population[0].fitness;
                }

                double gap = double(population[0].fitness - fitness_limit) / population[0].fitness * 100;
                if (gap < 0.01) {
                    break;
                }

                assert((int)population.size() == pop_size);

                fmax = INF;
                favg = 0;
                for (int i = 0; i < pop_size; i++) {
                    Individuo ind = population[i];
                    fmax = std::min(fmax, 1ll * ind.fitness);
                    favg += 1ll * ind.fitness;
                }
                favg /= pop_size;

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

                run_local_search(0, next_gen);
                run_crossover(elitism_size, next_gen);
                run_mutation(elitism_size, next_gen);


                std::sort(next_gen.begin(), next_gen.end());
                for (int i = elitism_size; i < pop_size; i++) {
                    population[i] = next_gen[i-elitism_size];
                }
                time(&at);
            }
        };

        population.assign(pop_size, Individuo(problem));
        int elitism_size = std::min(int(pop_size * elitism), 3);
        debug("creatining population");
        num_threads += 2;
        populate();
        num_threads -= 2;

        int gen = 0;
        std::sort(population.begin(), population.end());

        time_t start, nd, at;
        time(&start);
        time(&nd);
        time(&at);
        while ((double)(at - start) < main_tl) {
            double gap = double(population[0].fitness - fitness_limit) / population[0].fitness * 100;
            if (verbose) {
                std::cout << std::string(90, ' ') << '\r';
                std::cout.flush();
                std::cout << "Generation " << gen++ << ", Fitness =" << population[0].fitness << std::fixed << std::setprecision(6) << ", Gap: " << gap << "\n";
            }
            if (gap < 0.01) {
                break;
            }

            std::vector<Neighbor> vns_solutions(pop_size, Neighbor(problem));
            std::thread vns_thread, ga_thread;
            vns_thread = std::thread(run_vns, std::cref(population), std::cref(problem), std::ref(vns_solutions),
                                second_tl, fitness_limit, num_threads, pop_size);
            ga_thread = std::thread(run_ga, elitism_size);
            vns_thread.join();
            ga_thread.join();

            merge_solutions(vns_solutions);

            std::sort(population.begin(), population.end());
            time(&at);
        }

        if (verbose) {
            std::cout << std::string(90, ' ') << '\r';
            std::cout.flush();
            double gap = double(population[0].fitness - fitness_limit) / population[0].fitness * 100;
            std::cout << "Generation " << gen << ", Fitness =" << population[0].fitness << std::fixed << std::setprecision(6) << ", Gap: " << gap << "\n";
        }

        return {population[0], (double)(nd - start)};
    }

    void merge_solutions(std::vector<Neighbor> const &vns_solutions) {
        std::vector<Individuo> aux_vet(pop_size, Individuo(problem));
        for (int i = 0, l = 0, r = 0; i < pop_size; i++) {
            if (population[l].fitness < vns_solutions[i].fitness) {
                aux_vet[i] = population[l++];
            }
            else {
                aux_vet[i].dna = vns_solutions[r++].dna;
            }
            aux_vet[i].get_fitness();
        }
        std::swap(population, aux_vet);
    }
};