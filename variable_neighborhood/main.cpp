#include <iostream>
#include <assert.h>
#include <fstream>
#include "variable_neighborhood.cpp"
#include "neighbor.cpp"
#include "leasing.cpp"

// 1: my/path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 5);

    std::string dir = argv[1], filename = argv[2], version = argv[3];
    std::transform(version.begin(), version.end(), version.begin(), [](unsigned char c){ return std::toupper(c); });
    LeasingProblem problem(dir + filename, version);
    VNS model(problem);

    std::string fitness_limit = argv[4];

    auto [S, runtime] = model.run(60*15, std::min(std::max(problem.T * problem.V * problem.L * problem.K, 1000), 10000), true, std::stoll(fitness_limit));
    std::ofstream outfile;
    outfile.open("results.txt", std::ios_base::app);
    outfile << S << " time: " << runtime << '\n';
    outfile.close();
    std::cout << S << '\n' << runtime << '\n';



    return 0;
}
