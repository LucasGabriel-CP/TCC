#include <iostream>
#include <assert.h>
#include "variable_neighborhood.cpp"
#include "neighbor.cpp"
#include "leasing.cpp"

// 1: my/path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 4);

    std::string dir = argv[1], filename = argv[2], version = argv[3];
    std::transform(version.begin(), version.end(), version.begin(), [](unsigned char c){ return std::toupper(c); });
    LeasingProblem problem(dir + filename, version);
    VNS model(problem);

    auto [S, runtime] = model.run(60*10, problem.T * problem.V * problem.L * problem.K, true, 116);
    std::cout << S << '\n' << runtime << '\n';



    return 0;
}
