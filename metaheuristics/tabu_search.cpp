#include <vector>
#include <set>
#include <map>
#include <vector>
#include <queue>

#include "initial_solution.cpp"
#include "neighbor.cpp"
#include "problems/leasing.cpp"

struct TabuSearch {
    LeasingProblem problem;

    TabuSearch(LeasingProblem const &_problem) {
        problem = _problem;
    }

    std::pair<Neighbor, double> run(double time_limit = 3600, int r = 3, int k = 3, bool verbose = false, long long fitness_limit = -1) {
        time_t start, nd, at;
        std::set<long long> tabu_list;
        std::queue<std::pair<int, long long>> tabu_timer;

        // Solucao inicial
        time(&start);
        InitialSolution S0(problem);
        S0.build();
        time(&nd);
        Neighbor best_solution = Neighbor(problem, S0.dna, S0.fitness);
        Neighbor S = best_solution;
        bool stop_criteria = true;
        int cnt = 0;
        int save_k = k;
        for (int iter = 0; stop_criteria; iter++) { // roda ateh atingir o limite
            // Remove da lista tabu quem ja deu o tempo
            while ((int)tabu_list.size() > k) {
                tabu_list.erase(tabu_timer.front().second);
                tabu_timer.pop();
            }

            // Seleciona vizinhos
            std::vector<Neighbor> neighborhood(r, S);
            for (int i = 0; i < r; i++) {
                Neighbor &ng = neighborhood[i];
                ng.first_improve_ls();
            }

            // Pega o melhor vizinho que nao esta na lista
            std::sort(neighborhood.begin(), neighborhood.end());
            int id = 0;
            while (id < r && tabu_list.find(neighborhood[id].get_hashed_sol()) != tabu_list.end()) {
                id++;
            }
            if (id == r) {
                int x = rand_i() % r;
                if (neighborhood[x].fitness < S.fitness * 1.5) {
                    id = x;
                }
            }
            if (id != r) {
                tabu_list.insert(neighborhood[id].get_hashed_sol());
                tabu_timer.push({iter, neighborhood[id].get_hashed_sol()});
                S = neighborhood[id];
                if (neighborhood[id] < best_solution) {
                    time(&nd);
                    best_solution = neighborhood[id];
                    cnt = 0;
                    k = save_k;
                }
                else {
                    cnt++;
                }
            }
            if (cnt >= k) {
                if (rand() < 0.5) {
                    k = save_k;
                    tabu_list.clear();
                    tabu_timer = std::queue<std::pair<int, long long>>();
                }
                else {
                    k *= 1.5;
                }
                cnt = 0;
            }
            time(&at);
            double gap = double(best_solution.fitness - fitness_limit) / best_solution.fitness * 100;
            if (verbose) {
                std::cout << "\33[2k" << "Generation " << iter + 1 << ", Tabu list size: " << (int)tabu_list.size() << ", Cur Fitness: " << S.fitness
                    << ", Fitness = " << best_solution.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\r";
                std::cout.flush();
            }
            stop_criteria = double(at - start) < time_limit && gap > 0.005;
        }
        if (verbose) std::cout << '\n';
        if (best_solution.give_penalties()) best_solution.fix_solution();
        return {best_solution, double(nd - start)};
    }
};