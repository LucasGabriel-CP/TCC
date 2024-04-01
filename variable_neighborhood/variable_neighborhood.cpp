#pragma once

#include <iostream>
#include <vector>
#include <unordered_set>
#include "util/DebugTemplate.cpp"
#include "neighbor.cpp"
#include "leasing.cpp"
#include "initial_solution.cpp"

struct VNS {
    LeasingProblem problem;

    VNS(LeasingProblem const &_problem) {
        problem = _problem;
    }

    void sequential_change(Neighbor &S, Neighbor &SS, int &k, clock_t &nd) {
        if (SS < S) {
            nd = clock();
            std::swap(S, SS);
            k = 0;
        }
        else {
            k++;
        }
    }

    void pipe_change(Neighbor &S, Neighbor &SS, int &k, clock_t &nd) {
        if (SS < S) {
            nd = clock();
            std::swap(S, SS);
        }
        else {
            k++;
        }
    }

    void cyclic_change(Neighbor &S, Neighbor &SS, int &k, clock_t &nd) {
        if (SS < S) {
            nd = clock();
            std::swap(S, SS);
        }
        k++;
    }

    std::pair<Neighbor, double> run(double time_limit = 3600, int r = 3, bool verbose = false, long long fitness_limit = -1) {
        clock_t start = clock();
        InitialSolution S0(problem);
        S0.build();
        Neighbor S(problem, S0.dna, S0.fitness);
        clock_t nd;
        bool stop_criteria = bool(((double)(clock() - start) / CLOCKS_PER_SEC) < time_limit);
        int gen = 0;
        double gap = 0;
        std::unordered_set<long long> H;
        nd = (clock() - start) / CLOCKS_PER_SEC;
        while (stop_criteria) {
            Neighbor new_S = S;

            if (true) {
                gap = double(S.fitness - fitness_limit) / fitness_limit * 100;
                std::cout << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\n";
            }

            int k = 0;
            int cnt = 0;
            if ((int)H.size() >= 10000000) {
                H.clear();
            }
            while(k < r && stop_criteria) {
                Neighbor SS = new_S;
                long long hash = SS.get_hashed_sol();
                int tried = 0;
                SS.shake();
                while (H.find(hash) != H.end() && tried < 100) {
                    SS = new_S;
                    SS.shake();
                    hash = SS.get_hashed_sol();
                    tried++;
                }
                if (tried != 100) cnt++;
                H.insert(hash);
                // if (SS.shake_mutate()) cnt++;
                SS.local_searchv2();
                
                // Select change
                sequential_change(new_S, SS, k, nd);

                if (verbose) {
                    gap = double(SS.fitness - fitness_limit) / fitness_limit * 100;
                    std::cout << "\33[2k" << "Neighbor " << k << ", Fitness = " << SS.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\r";
                    std::cout.flush();
                }
                stop_criteria = bool(((double)(clock() - start) / CLOCKS_PER_SEC) < time_limit) && gap > 0.005;
            }
            if (verbose) std::cout << '\n';
            debug(cnt);
            if (new_S < S) {
                std::swap(new_S, S);
            }
            stop_criteria = bool(((double)(clock() - start) / CLOCKS_PER_SEC) < time_limit) && gap > 0.005;
            gen++;
        }

        assert(S.give_penalties() == 0);

        if (true) {
            gap = double(S.fitness - fitness_limit) / fitness_limit * 100;
            std::cout << std::string(100, ' ') << '\r';
            std::cout.flush();
            std::cout << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\n";
            std::cout.flush();
        }
        return {S, double(nd - start) / CLOCKS_PER_SEC};
    }
};