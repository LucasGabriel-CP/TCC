#pragma once

#include <iostream>
#include <vector>
#include <unordered_set>
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

    Neighbor vnd(Neighbor &S0, int r = 3, long long fitness_limit = -1) {
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
            SS.first_improve_ls();
            
            // Select change
            sequential_change(new_S, SS, k, nd);

            time(&at);
            gap = double(SS.fitness - fitness_limit) / SS.fitness * 100;
            stop_criteria = gap > 0.005;
        }
        if (new_S < S) {
            std::swap(new_S, S);
        }
        return S;
    }

    std::pair<Neighbor, double> run(double time_limit = 3600, int r = 3, bool verbose = false, long long fitness_limit = -1, int lmax = 3) {
        time_t start, nd, at;
        time(&start);
        InitialSolution S0(problem);
        S0.build();
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
                SS.shake();
                SS = vnd(SS, lmax, fitness_limit);
                
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

        // if (S.give_penalties()) {
        //     S.fix_solution();
        // }

        if (verbose) {
            gap = double(S.fitness - fitness_limit) / S.fitness * 100;
            std::cout << std::string(100, ' ') << '\n';
            std::cout << "Generation " << gen << ", Fitness = " << S.fitness << ", Gap:" << std::fixed << std::setprecision(4) << gap << "\n";
            std::cout.flush();
        }
        return {S, double(nd - start)};
    }
};