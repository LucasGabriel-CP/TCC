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
    Individuo generate_individuo(std::vector<int> const &order) {
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
    void populate_partially(int l, int r, std::vector<int> const &order) {
        r = std::min(r, (int)population.size()-1);
        for (int i = l; i <= r; i++) {
            population[i] = generate_individuo(order);
        }
    }

    /**
     * Function to initialize the process of population.
     * @param num_threads Number of threads to be created.
    */
    void populate(int num_threads) {
        // Problem specific stuff, this shall change.
        std::vector<int> order(num_t);
        // std::iota(order.begin(), order.end(), 0);
        // std::sort(order.begin(), order.end(), [&](const int &a, const int &b) {
        //     double x1 = turbinies[a].pos.x - substation.x;
        //     double y1 = turbinies[a].pos.y - substation.y;
        //     double x2 = turbinies[b].pos.x - substation.x;
        //     double y2 = turbinies[b].pos.y - substation.y;
        //     return std::atan2(x1, y1) < std::atan2(x2, y2);
        // });

        std::vector<std::thread> threads(num_threads);
        int l = 0, step = (int)population.size()/num_threads;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads-1) step = (int)population.size(); // Give remainder to last thread, prob better to give this to the first one.
            int r = l + step-1;
            threads[i] = std::thread(&Evolution::populate_partially, this, l, r, std::cref(order));
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
        // parents.erase(std::unique(parents.begin(), parents.end()), parents.end());

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
     * Function to handle the evolutionary process.
    */
    std::pair<Individuo, double> run_evo(
        double fitness_limit = -1, double elitism = 0.1, int pop_size = 100,
        double time_limit = 3600, bool verbose = false, int num_threads = 4
    ) {
        population.resize(pop_size);
        int elitism_size = std::min(int(pop_size * elitism), 3);
        int start = clock();
        populate(num_threads);

        double best = DOUBLE_INF;
        int gen = 0;
        while ((double)(clock() - start) / CLOCKS_PER_SEC < time_limit) {
            std::sort(population.begin(), population.end());
            best = std::min(best, population[0].fitness);

            if (fabs(best - fitness_limit) < EPS) {
                break;
            }

            while ((int)population.size() > pop_size) population.pop_back();
            assert((int)population.size() == pop_size);
            double worst = population.back().fitness;

            double fmax = DOUBLE_INF, favg = 0;
            for (Individuo &ind: population) {
                fmax = std::min(fmax, ind.fitness);
                favg += ind.fitness;
            }
            favg /= pop_size;

            if (verbose) {
                double gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2) * 100;
                std::cout << std::string(80, ' ') << '\r';
                std::cout.flush();
                std::cout << "Generation " << gen << ", Fitness: min="
                    << std::fixed << std::setprecision(6) << fmax << ", mean=" << favg
                    << ", max=" << worst << ", Gap:" << gap << "\r";
                std::cout.flush();
            }
            std::vector<int> weigths(pop_size);
            int cur_weight = pop_size + 5;
            double ant = best;
            for (int i = 0; i < pop_size; i++) {
                if (f_cmp(population[i].fitness, ant) > 0) {
                    ant = population[i].fitness;
                    cur_weight--;
                }
                weigths[i] = cur_weight;
            }
            d = std::discrete_distribution<int>(weigths.begin(), weigths.end());

            std::vector<std::thread> threads(num_threads);
            std::vector<Individuo> next_gen(pop_size);
            int l = elitism_size, step = (pop_size - elitism_size) / num_threads;
            for (int i = 0; i < num_threads; i++) {
                int r;
                if (i == num_threads - 1) r = pop_size - 1;
                else r = l + step - 1;
                threads[i] = std::thread(&Evolution::run_crossover_partially, this, fmax, favg, l, r, std::ref(next_gen));
                l += step;
            }
            for (std::thread &th: threads) th.join();
            
            l = elitism_size, step = (pop_size - elitism_size) / num_threads;
            for (int i = 0; i < num_threads; i++) {
                int r;
                if (i == num_threads - 1) r = pop_size - 1;
                else r = l + step - 1;
                threads[i] = std::thread(&Evolution::run_mutation_partially, this, fmax, favg, l, r, std::ref(next_gen));
                l += step;
            }
            for (std::thread &th: threads) th.join();

            for (int i = 0; i < elitism_size; i++) next_gen[i] = population[i];
            population = std::move(next_gen);

            gen++;
        }


        if (verbose) {
            double gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2) * 100;
            std::cout << "Generation " << gen << ", Fitness =" << std::fixed << std::setprecision(6) << best << ", Gap: " << gap << "\n";
        }

        int nd = (clock() - start) / CLOCKS_PER_SEC;

        return {population[0], nd};
    }
};
