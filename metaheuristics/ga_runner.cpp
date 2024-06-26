#include <iostream>
#include "genetic_algorithm.cpp"
#include "problems/leasing.cpp"

// 1: path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    debug(argv);
    assert(argc == 5);

    std::string dir = argv[1], version = argv[3];
    LeasingProblem problem(dir + argv[2], version);


    std::string fitness_limit = argv[4];
    double elitism = .1;
    int pop_size = 32;
    int num_threads = 16;
    double time_limit = 60*15;
    bool verbose = true;
    Evolution ga(problem, pop_size, num_threads);

    auto [ind, t] = ga.run_evo(std::stoll(fitness_limit), elitism, time_limit, verbose);
    std::ofstream outfile;
    time_t curr_time = std::time(NULL);
    std::string curr_date = std::ctime(&curr_time);
    curr_date.pop_back();
    outfile.open("logs/ga.log", std::ios_base::app);
    outfile << curr_date << ' ' << argv[2] << ' ' << ind << " time: " << t << '\n';
    outfile.close();
    std::cout << ind << '\n';
    std::cout << std::fixed << std::setprecision(4) << t << '\n';


    return 0;
}
