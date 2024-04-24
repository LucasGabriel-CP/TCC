#include <iostream>
#include <assert.h>
#include <fstream>
#include <thread>
#include "genetic_algorithm.cpp"
#include "problems/leasing.cpp"
#include "variable_neighborhood.cpp"
#include "neighbor.cpp"

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


// 1: path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    debug(argv);
    assert(argc == 5);

    std::string dir = argv[1], version = argv[3];
    LeasingProblem problem(dir + argv[2], version);


    std::string fitness_limit = argv[4];
    double elitism = .1;
    int pop_size = 32;
    int num_threads = 14;
    bool verbose = true;
    
    double main_tl = 60 * 15;
    double second_tl = 60 * 3;

    Evolution ga(problem, pop_size, num_threads);

    auto [ind, t] = ga.run_hibrid(std::stoll(fitness_limit), elitism, main_tl, second_tl, verbose);
    std::ofstream outfile;
    time_t curr_time = std::time(NULL);
    std::string curr_date = std::ctime(&curr_time);
    curr_date.pop_back();
    outfile.open("logs/hibrid.log", std::ios_base::app);
    outfile << curr_date << ' ' << argv[2] << ' ' << ind << " time: " << t << '\n';
    outfile.close();
    std::cout << ind << '\n';
    std::cout << std::fixed << std::setprecision(4) << t << '\n';

    return 0;
}
