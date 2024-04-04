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
        for (int iter = 0; stop_criteria; iter++) { // roda ateh atingir o limite
            // Remove da lista tabu quem ja deu o tempo
            if ((int)tabu_list.size() == k) {
                tabu_list.erase(tabu_timer.front().second);
                tabu_timer.pop();
            }

            // Seleciona vizinhos
            std::vector<Neighbor> neighborhood(r, S);
            for (Neighbor &ng: neighborhood) {
                // ng.shake();
                ng.best_improve_ls();
            }

            // Pega o melhor vizinho que nao esta na lista
            std::sort(neighborhood.begin(), neighborhood.end());
            int id = 0;
            while (id < r && tabu_list.find(neighborhood[id].get_hashed_sol()) != tabu_list.end()) {
                id++;
            }
            if (id != r) {
                tabu_list.insert(neighborhood[id].get_hashed_sol());
                tabu_timer.push({iter, neighborhood[id].get_hashed_sol()});
                S = neighborhood[id];
                if (neighborhood[id] < best_solution) {
                    time(&nd);
                    best_solution = neighborhood[id];
                }
            }
            time(&at);
            double gap = double(best_solution.fitness - fitness_limit) / best_solution.fitness * 100;
            if (verbose) {
                std::cout << "\33[2k" << "Generation " << iter << ", Tabu list size: " << (int)tabu_list.size() << ", Cur Fitness: " << S.fitness
                    << ", Fitness = " << best_solution.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\r";
                std::cout.flush();
            }
            stop_criteria = double(at - start) < time_limit && gap > 0.005;
        }
        if (verbose) std::cout << '\n';
        return {best_solution, double(nd - start)};
    }
};