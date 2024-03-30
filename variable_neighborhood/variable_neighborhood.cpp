#pragma once

#include <iostream>
#include <vector>
#include "util/DebugTemplate.cpp"
#include "neighbor.cpp"
#include "leasing.cpp"
#include "initial_solution.cpp"

struct VNS {
    LeasingProblem problem;

    VNS(LeasingProblem const &_problem) {
        problem = _problem;
    }

    void sequential_change(Neighbor &S, Neighbor &SS, int &k) {
        if (SS < S) {
            std::swap(S, SS);
            k = 0;
        }
        else {
            k++;
        }
    }

    void pipe_change(Neighbor &S, Neighbor &SS, int &k) {
        if (SS < S) {
            std::swap(S, SS);
        }
        else {
            k++;
        }
    }

    void cyclic_change(Neighbor &S, Neighbor &SS, int &k) {
        if (SS < S) {
            std::swap(S, SS);
        }
        k++;
    }

    std::pair<Neighbor, double> run(double time_limit = 3600, int r = 3, bool verbose = false, long long fitness_limit = -1) {
        double start = clock();
        InitialSolution S0(problem);
        S0.build();
        Neighbor S(problem, S0.dna, S0.fitness);
        double nd;
        bool stop_criteria = bool((double)(clock() - start) / CLOCKS_PER_SEC < time_limit);
        int gen = 0;
        while (stop_criteria) {
            Neighbor new_S = S;
            nd = (clock() - start) / CLOCKS_PER_SEC;

            if (verbose) {
                double gap = (S.fitness - fitness_limit) / ((S.fitness + fitness_limit) / 2.0) * 100;
                std::cout << std::string(100, ' ') << '\r';
                std::cout.flush();
                std::cout << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(6) << gap << "\r";
                std::cout.flush();
            }

            new_S.shake();
            int k = 0;
            while(k < r) {
                Neighbor SS = new_S;
                SS.local_search();
                
                // Select change
                pipe_change(new_S, SS, k);
            }
            if (new_S < S) {
                std::swap(new_S, S);
                stop_criteria = bool((double)(clock() - start) / CLOCKS_PER_SEC < time_limit);
            }
            else {
                stop_criteria = false;
            }
            gen++;
        }

        if (verbose) {
            double gap = (S.fitness - fitness_limit) / ((S.fitness + fitness_limit) / 2.0) * 100;
            std::cout << std::string(100, ' ') << '\r';
            std::cout.flush();
            std::cout << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(6) << gap << "\n";
            std::cout.flush();
        }
        return {S, nd};
    }
};