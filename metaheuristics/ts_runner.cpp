#include <iostream>
#include <assert.h>
#include <fstream>
#include <thread>
#include "tabu_search.cpp"
#include "neighbor.cpp"
#include "problems/leasing.cpp"

struct params {
    int r, k;
    long long fitness;
    params() { };
    params(int _r, int _k, long long _fitness) {
        r = _r; _k = k;
        fitness = _fitness;
    }
    friend bool operator < (params const &lhs, params const &rhs){
        return lhs.fitness < rhs.fitness;
    }

    friend std::ostream &operator <<(std::ostream &os, const params &sol) {
        os << "Fitness: " << sol.fitness << ", r: " << sol.r << ", k: " << sol.k;
        return os;
    }
};

// 1: my/path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 5);

    std::string dir = argv[1], filename = argv[2], version = argv[3];
    std::transform(version.begin(), version.end(), version.begin(), [](unsigned char c){ return std::toupper(c); });
    LeasingProblem problem(dir + filename, version);
    std::string fitness_limit = argv[4];


    // TabuSearch model(problem);
    // auto [S, runtime] = model.run(60*15, 50, 75, true, std::stoll(fitness_limit));
    // std::ofstream outfile;
    // outfile.open("./logs/ts.log", std::ios_base::app);
    // outfile << filename << ' ' << S << " time: " << runtime << '\n';
    // outfile.close();
    // std::cout << S << '\n' << runtime << '\n';

    params best_solution(INF, -1, -1);
    int sr = 1, sk = 1;
    while (sr <= 50) {
        int solution_number = 16;
        std::vector<std::thread> threads(solution_number);
        std::vector<params> all_fitness(solution_number);
        auto vns_runner = [&](int i, int r, int k) {
            TabuSearch model(problem);
            auto [S, runtime] = model.run(60*30, r, k, false, std::stoll(fitness_limit));
            mtx.lock();
                if (S.give_penalties() != 0) std::cout << "Invalid solution\n";
                std::cout << S << '\n' << runtime << '\n';
                all_fitness[i] = {r, k, S.fitness};
            mtx.unlock();
        };
        for (int i = 0; i < solution_number; i++) {
            threads[i] = std::thread(vns_runner, i, sr, sk);
            sk = (sk + 1) % 100;
            if (!sk) {
                sk++;
                sr++;
            }
        }

        for (std::thread &th: threads) th.join();

        for (auto f: all_fitness) {
            if (f < best_solution) {
                best_solution = f;
            }
        }
    }
    std::cout << "Best fitness overall: " << best_solution << '\n';
    std::ofstream outfile;
    outfile.open("logs/ts.log", std::ios_base::app);
    outfile << argv[2] << ' ' << best_solution << '\n';
    outfile.close();

    return 0;
}
