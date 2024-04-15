#include <iostream>
#include <assert.h>
#include <fstream>
#include <thread>
#include "tabu_search.cpp"
#include "neighbor.cpp"
#include "problems/leasing.cpp"

struct params {
    double runtime;
    long long fitness;
    params() { }
    params(double _runtime, long long _fitness) {
        runtime = _runtime; fitness = _fitness;
    }
    friend bool operator < (params const &lhs, params const &rhs){
        if (lhs.fitness != rhs.fitness) return lhs.fitness < rhs.fitness;
        return lhs.runtime < rhs.runtime;
    }

    friend std::ostream &operator <<(std::ostream &os, const params &sol) {
        os << "Fitness: " << sol.fitness << ", runtime: " << sol.runtime;
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

    int sr = 15, sk = 48;
    int solution_number = 16;
    std::vector<std::thread> threads(solution_number);
    std::vector<params> all_fitness(solution_number);
    auto vns_runner = [&](int i, int r, int k) {
        TabuSearch model(problem);
        auto [S, runtime] = model.run(60*15, r, k, false, std::stoll(fitness_limit));
        mtx.lock();
            if (S.give_penalties() != 0) std::cout << "Invalid solution\n";
            all_fitness[i] = {runtime, S.fitness};
            std::cout << all_fitness[i] << '\n' << runtime << '\n';
        mtx.unlock();
    };
    for (int i = 0; i < solution_number; i++) {
        threads[i] = std::thread(vns_runner, i, sr, sk);
    }

    for (std::thread &th: threads) th.join();

    params best_solution = *std::min_element(all_fitness.begin(), all_fitness.end());
    time_t curr_time = std::time(NULL);
    std::string curr_date = std::ctime(&curr_time);
    curr_date.pop_back();
    std::cout << "Best fitness overall: " << best_solution << '\n';
    std::ofstream outfile;
    outfile.open("logs/ts.log", std::ios_base::app);
    outfile << curr_date << ' ' << argv[2] << ' ' << best_solution << '\n';
    outfile.close();

    return 0;
}
