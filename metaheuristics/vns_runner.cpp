#include <iostream>
#include <assert.h>
#include <fstream>
#include <thread>
#include "variable_neighborhood.cpp"
#include "neighbor.cpp"
#include "problems/leasing.cpp"

// 1: my/path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 5);

    std::string dir = argv[1], filename = argv[2], version = argv[3];
    std::transform(version.begin(), version.end(), version.begin(), [](unsigned char c){ return std::toupper(c); });
    LeasingProblem problem(dir + filename, version);
    std::string fitness_limit = argv[4];


    // VNS model(problem);
    // auto [S, runtime] = model.run(60*15, 50, true, std::stoll(fitness_limit));
    // std::ofstream outfile;
    // outfile.open("results.txt", std::ios_base::app);
    // outfile << S << " time: " << runtime << '\n';
    // outfile.close();
    // std::cout << S << '\n' << runtime << '\n';

    return 0;

    int solution_number = 12;
    std::vector<std::thread> threads(solution_number);
    std::vector<long long> all_fitness(solution_number);
    auto vns_runner = [&](int i) {
        VNS model(problem);
        auto [S, runtime] = model.run(60*15, 50, false, std::stoll(fitness_limit));
        mtx.lock();
            if (S.give_penalties() != 0) std::cout << "Invalid solution\n";
            std::cout << S << '\n' << runtime << '\n';
            all_fitness[i] = S.fitness;
        mtx.unlock();
    };
    for (int i = 0; i < solution_number; i++) {
        threads[i] = std::thread(vns_runner, i);
    }

    for (std::thread &th: threads) th.join();

    long long best = INF;
    for (auto f: all_fitness) best = std::min(best, f);
    std::cout << "Best fitness overall: " << best << '\n';
    std::ofstream outfile;
    outfile.open("results.txt", std::ios_base::app);
    outfile << argv[2] << ' ' << best << '\n';
    outfile.close();

    return 0;
}
