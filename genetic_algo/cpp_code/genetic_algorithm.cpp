#include "util/Individuo.cpp"
#include "util/Constants.cpp"

struct Evolution {
    std::vector<Individuo> population;
    std::vector<Individuo> next_gen;
    std::discrete_distribution<int> d;
    Evolution() { }

    Individuo generate_individuo() {
        Individuo new_ind = Individuo();

        return new_ind;
    }

    void populate(int l, int r) {
        r = std::min(r, (int)population.size());
        for (int i = l; i <= r; i++) {
            population[i] = generate_individuo();
        }
    }

    void populate_partially(int l, int r) {
        r = std::min(r, (int)population.size());
        for (int i = l; i <= r; i++) {
            population[i] = generate_individuo();
        }
    }

    void populate() {
        std::vector<std::thread> threads(5);
        int l = 0, step = (int)population.size()/5;
        for (int i = 0; i < 5; i++) {
            if (i == 4) step = (int)population.size();
            int r = l + step-1;
            threads[i] = std::thread(&populate_partially, this, l, r);
            l += step;
        }
        for (std::thread &th: threads) th.join();
    }

    std::pair<Individuo, Individuo> select_parents() {
        std::vector<int> parents;
        for (int i = 0; i < 5; i++) parents.push_back(d(g));
        std::sort(parents.begin(), parents.end(), [&](const int &a, const int &b){
            return population[a] < population[b];
        });
        return {population[parents[0]], Individuo(population[parents[1]])};
    }

    float crossover_prob(float fmax, float favg, float fat) {
        float k3 = 0.5;
        float k1 = rand();
        if (fat <= favg) {
            return k3;
        }
        return k1 * (favg - fmax) / std::max(fat - fmax, 1.0F);
    }

    float mutation_prob(float fmax, float favg, float fat) {
        float k4 = 0.5;
        float k2 = rand();
        if (fat <= favg) {
            return k4;
        }
        return k2 * (favg - fmax) / std::max(fat - fmax, 1.0F);
    }

    void run_evo_partially(float fmax, float favg, int l, int r) {
        for (int i = l; i < r; i++) {
            std::pair<Individuo, Individuo> parents = select_parents();
            if (rand() < crossover_prob(fmax, favg, population[i].fitness)) {
                auto [child_1, child_2] = Individuo::cross(parents.first, parents.second);
                mtx.lock();
                    next_gen.push_back(child_1);
                    next_gen.push_back(child_2);
                mtx.unlock();
            }
            if (rand() < mutation_prob(fmax, favg, population[i].fitness)) {
                population[i].mutate();
                population[i].get_fitness();
            }
            next_gen[i] = population[i];
        }
    }

    std::pair<Individuo, float> run_evo(
        float fitness_limit = -1, float elitism = 0.1, int pop_size = 100,
        float time_limit = 3600, bool verbose = false, int num_threads = 4
    ) {
        population.resize(pop_size);
        int elitism_size = std::min(int(pop_size * elitism), 3);
        int start = clock();
        populate();
        float best = FLOAT_INF;
        int gen = 0;
        while ((double)(clock() - start) / CLOCKS_PER_SEC < time_limit) {
            std::sort(population.begin(), population.end());
            best = std::min(best, population[0].fitness);


            if (fabs(best - fitness_limit) < EPS) {
                break;
            }

            while ((int)population.size() > pop_size) population.pop_back();
            float worst = population.back().fitness;

            float fmax = FLOAT_INF, favg = 0;
            for (Individuo &ind: population) {
                fmax = std::min(fmax, ind.fitness);
                favg += ind.fitness;
            }
            favg /= pop_size;

            if (verbose) {
                float gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2) * 100;
                std::cout << "Generation " << gen << ", Fitness: min="
                    << std::fixed << std::setprecision(6) << fmax << ", mean=" << favg
                    << ", max=" << worst << ", Gap:" << gap << "\r";
                std::cout.flush();
            }

            std::vector<int> weigths(pop_size);
            for (int i = 0; i < pop_size; i++) {
                weigths[i] = 1 + worst - population[i].fitness;
            }
            d = std::discrete_distribution<int>(weigths.begin(), weigths.end());

            next_gen = std::vector<Individuo>(pop_size);

            for (int i = 0; i < elitism_size; i++) next_gen[i] = population[i];

            std::vector<std::thread> threads(num_threads);
            int l = 0, step = pop_size / num_threads;
            for (int i = 0; i < num_threads; i++) {
                if (i == num_threads - 1) step = pop_size;
                int r = l + step - 1;
                threads[i] = std::thread(&run_evo_partially, this, fmax, favg, l, r);
                l += step;
            }

            for (std::thread &th: threads) th.join();

            std::swap(population, next_gen);
            gen++;
        }

        if (verbose) {
            float gap = (population[0].fitness - fitness_limit) / ((population[0].fitness + fitness_limit) / 2) * 100;
            std::cout << "Generation " << gen << ", Fitness =" << std::fixed << std::setprecision(6) << best << ", Gap: " << gap << "\n";
        }

        int nd = (clock() - start) / CLOCKS_PER_SEC;

        return {population[0], nd};
    }
};
