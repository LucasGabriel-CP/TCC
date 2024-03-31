#include <iostream>
#include "genetic_algorithm.cpp"
#include "leasing.cpp"

// 1: path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    debug(argv);
    assert(argc == 4);

    std::string dir = argv[1], version = argv[3];
    LeasingProblem problem(dir + argv[2], version);


    int fitness_limit = 23168;
    double elitism = .1;
    int pop_size = 200;
    int num_threads = 10;
    double time_limit = 60*15;
    bool verbose = true;
    Evolution ga(problem, pop_size, num_threads);

    auto [ind, t] = ga.run_evo(fitness_limit, elitism, time_limit, verbose);
    std::cout << ind << '\n';
    std::cout << std::fixed << std::setprecision(5) << t << '\n';


    return 0;
}
