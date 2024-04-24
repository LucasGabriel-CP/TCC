#pragma once

#include <iostream>
#include <vector>
#include <unordered_set>
#include <queue>
#include "util/DebugTemplate.cpp"
#include "neighbor.cpp"
#include "problems/leasing.cpp"
#include "initial_solution.cpp"

struct VNS {
    LeasingProblem problem;

    VNS(LeasingProblem const &_problem) {
        problem = _problem;
    }

    void sequential_change(Neighbor &S, Neighbor &SS, int &k, time_t &nd) {
        if (SS < S) {
            time(&nd);
            std::swap(S, SS);
            k = 0;
        }
        else {
            k++;
        }
    }

    void pipe_change(Neighbor &S, Neighbor &SS, int &k, time_t &nd) {
        if (SS < S) {
            time(&nd);
            std::swap(S, SS);
        }
        else {
            k++;
        }
    }

    void cyclic_change(Neighbor &S, Neighbor &SS, int &k, time_t &nd) {
        if (SS < S) {
            time(&nd);
            std::swap(S, SS);
        }
        k++;
    }

    Neighbor vnd(Neighbor &S0, int r = 3, long long fitness_limit = -1, double time_limit = 3600) {
        time_t start, nd, at;
        time(&start);
        Neighbor S(problem, S0.dna, S0.fitness);
        double gap = 0;
        time(&nd);
        bool stop_criteria = true;
        int k = 0;
        Neighbor new_S = S;
        while(k < r && stop_criteria) {
            Neighbor SS = new_S;
            // SS.get_neighbor(k+1);
            SS.first_improve_ls(k+1);

            // Select change
            sequential_change(new_S, SS, k, nd);

            time(&at);
            gap = double(SS.fitness - fitness_limit) / SS.fitness * 100;
            stop_criteria = bool(double(at - start) < time_limit) && gap > 0.01;
        }
        if (new_S < S) {
            std::swap(new_S, S);
        }
        return S;
    }

    std::pair<Neighbor, double> run(
        double time_limit = 3600, int r = 3, bool verbose = false,
        long long fitness_limit = -1, int lmax = 3, bool use_tabu = false, int tabu_k = 3
    ) {
        time_t start, nd, at;
        int save_k = tabu_k;
        std::set<long long> tabu_list;
        std::queue<long long> tabu_timer;
        int gen = 0;
        double gap = 0;

        time(&start);
        InitialSolution S0(problem);
        S0.random_build();
        Neighbor S(problem, S0.dna, S0.fitness);
        time(&nd);
        time(&at);
        bool stop_criteria = bool(double(nd - start) < time_limit);
        int cnt = 0;
        debug(use_tabu, tabu_k);
        while (stop_criteria) {
            Neighbor new_S = S;

            if (verbose) {
                std::cout << "\33[2k" << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\n";
                std::cout.flush();
            }

            int k = 0;
            while(k < r && stop_criteria) {
                Neighbor SS = new_S;
                int _i = 0;
                if (use_tabu) {
                    while ((int)tabu_list.size() > k) {
                        tabu_list.erase(tabu_timer.front());
                        tabu_timer.pop();
                    }
                    std::vector<Neighbor> neighborhood(k+1, SS);
                    for (int i = 0; i < k+1; i++) {
                        neighborhood[i].shake(k+1);
                    }

                    std::sort(neighborhood.begin(), neighborhood.end());
                    int id = 0;
                    while (id < (k + 1) && tabu_list.find(neighborhood[id].get_hashed_sol()) != tabu_list.end()) {
                        id++;
                    }
                    if (id == (k + 1)) {
                        int x = rand_i() % (k + 1);
                        if (neighborhood[x].fitness < S.fitness * 1.5) {
                            id = x;
                        }
                    }
                    if (id != (k + 1)) {
                        SS = neighborhood[id];
                        long long hash_val = SS.get_hashed_sol();
                        tabu_list.insert(hash_val);
                        tabu_timer.push(hash_val);
                    }
                }
                else {
                    SS.shake(k+1);
                }
                SS = vnd(SS, lmax, fitness_limit, time_limit - double(at - start));
                // SS.first_improve_ls(k+1);
                
                // Select change
                sequential_change(new_S, SS, k, nd);

                if (k) {
                    cnt++;
                }
                else {
                    cnt = 0;
                    tabu_k = save_k;
                }

                time(&at);
                gap = double(SS.fitness - fitness_limit) / fitness_limit * 100;
                stop_criteria = bool(double(at - start) < time_limit) && gap > 0.01;

                if (cnt >= tabu_k) {
                    if (rand() < 0.5) {
                        tabu_k = save_k;
                        tabu_list.clear();
                        tabu_timer = std::queue<long long>();
                    }
                    else {
                        tabu_k *= 1.5;
                    }
                    cnt = 0;
                }
                if (verbose) {
                    std::cout << "\33[2k" << "NeiN " << _i << ", Fitness = " << SS.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\n";
                    std::cout.flush();
                }
            }
            if (new_S < S) {
                std::swap(new_S, S);
            }
            gap = double(S.fitness - fitness_limit) / S.fitness * 100;
            stop_criteria = bool(double(at - start) < time_limit) && gap > 0.005;
            gen++;
        }

        if (S.give_penalties()) {
            S.fix_solution();
        }

        if (verbose) {
            gap = double(S.fitness - fitness_limit) / S.fitness * 100;
            std::cout << std::string(100, ' ') << '\n';
            std::cout << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\n";
            std::cout.flush();
        }
        return {S, double(nd - start)};
    }

    std::pair<Neighbor, double> run(Neighbor const &S0, double time_limit = 3600, int r = 3, bool verbose = false, long long fitness_limit = -1, int lmax = 3) {
        time_t start, nd, at;
        time(&start);
        Neighbor S(problem, S0.dna, S0.fitness);
        int gen = 0;
        double gap = 0;
        time(&nd);
        time(&at);
        bool stop_criteria = bool(double(nd - start) < time_limit);
        while (stop_criteria) {
            Neighbor new_S = S;

            if (verbose) {
                std::cout << "\33[2k" << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\r";
                std::cout.flush();
            }

            int k = 0;
            while(k < r && stop_criteria) {
                Neighbor SS = new_S;
                SS.shake(k+1);
                // SS = vnd(SS, lmax, fitness_limit);
                SS.first_improve_ls(k+1);
                
                // Select change
                sequential_change(new_S, SS, k, nd);

                time(&at);
                gap = double(SS.fitness - fitness_limit) / fitness_limit * 100;
                stop_criteria = bool(double(at - start) < time_limit) && gap > 0.005;
            }
            if (new_S < S) {
                std::swap(new_S, S);
            }
            gap = double(S.fitness - fitness_limit) / S.fitness * 100;
            stop_criteria = bool(double(at - start) < time_limit) && gap > 0.005;
            gen++;
        }

        if (S.give_penalties()) {
            S.fix_solution();
        }

        if (verbose) {
            gap = double(S.fitness - fitness_limit) / S.fitness * 100;
            std::cout << std::string(100, ' ') << '\n';
            std::cout << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\n";
            std::cout.flush();
        }
        return {S, double(nd - start)};
    }

    void run_multisol(std::vector<Neighbor> &solutions, int l_idx, int r_idx, double time_limit = 3600, int r = 3, long long fitness_limit = -1) {
        Neighbor S(problem, solutions[l_idx].dna, solutions[l_idx].fitness);
        double gap = 0;
        time_t start, nd, at;
        time(&start);
        time(&nd);
        time(&at);
        bool stop_criteria = bool(double(nd - start) < time_limit);
        while (stop_criteria) {
            Neighbor new_S = S;

            int k = 0;
            while(k < r && stop_criteria) {
                Neighbor SS = new_S;
                SS.shake(k+1);
                SS = vnd(SS, 48, fitness_limit, time_limit - double(at - start));
                // SS = vnd(SS, 10, fitness_limit, time_limit - double(at - start));
                // SS.first_improve_ls();
                
                // Select change
                sequential_change(new_S, SS, k, nd);

                time(&at);
                gap = double(SS.fitness - fitness_limit) / fitness_limit * 100;
                stop_criteria = bool(double(at - start) < time_limit) && gap > 0.009;
            }
            int id = std::upper_bound(solutions.begin() + l_idx, solutions.begin() + r_idx, new_S) - solutions.begin();
            if (id <= r_idx) {
                solutions[id] = new_S;
            }
            if (new_S < S) {
                std::swap(new_S, S);
            }
            gap = double(S.fitness - fitness_limit) / S.fitness * 100;
            stop_criteria = bool(double(at - start) < time_limit) && gap > 0.009;
        }
    }
};