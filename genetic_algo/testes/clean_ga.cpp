#include <iostream>
#include <thread>
#include <assert.h>
#include <iomanip>
#include <chrono>
#include "Individuo.cpp"
#include "util/Constants.cpp"

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
    void populate_partially(int l, int r, std::vector<Individuo> &initial_solution) {
        r = std::min(r, (int)initial_solution.size()-1);
        for (int i = l; i <= r; i++) {
            initial_solution[i] = generate_individuo();
        }
    }

    /**
     * Function to initialize the process of population.
     * @param num_threads Number of threads to be created.
    */
    void populate(int num_threads) {
        std::vector<std::thread> threads(num_threads);
        std::vector<Individuo> initial_solution((int)population.size());
        int l = 0, step = (int)population.size()/num_threads;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads-1) step = (int)population.size(); // Give remainder to last thread, prob better to give this to the first one.
            int r = l + step-1;
            threads[i] = std::thread(&Evolution::populate_partially, this, l, r, std::ref(initial_solution));
            l += step;
        }
        for (std::thread &th: threads) th.join();
        for (int i = 0; i < (int)population.size(); i++) {
            population[i] = initial_solution[i];
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
