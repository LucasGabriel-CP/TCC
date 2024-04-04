#include <iostream>
#include <assert.h>
#include <fstream>
#include <thread>
#include "tabu_search.cpp"
#include "neighbor.cpp"
#include "problems/leasing.cpp"

// 1: my/path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 5);

    std::string dir = argv[1], filename = argv[2], version = argv[3];
    std::transform(version.begin(), version.end(), version.begin(), [](unsigned char c){ return std::toupper(c); });
    LeasingProblem problem(dir + filename, version);
    std::string fitness_limit = argv[4];


    TabuSearch model(problem);
    auto [S, runtime] = model.run(60*15, 50, 25, true, std::stoll(fitness_limit));
    std::ofstream outfile;
    outfile.open("ts.log", std::ios_base::app);
    outfile << filename << ' ' << S << " time: " << runtime << '\n';
    outfile.close();
    std::cout << S << '\n' << runtime << '\n';

    return 0;
}
