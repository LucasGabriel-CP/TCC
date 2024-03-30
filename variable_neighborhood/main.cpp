#include <iostream>
#include <assert.h>
#include "variable_neighborhood.cpp"
#include "neighbor.cpp"
#include "leasing.cpp"

// 1: my/path/, 2: instance_filename.txt
int main(int argc, char *argv[]) {
    assert(argc == 3);

    std::string dir = argv[1], filename = argv[2];
    LeasingProblem problem(dir + filename);
    VNS model(problem);

    auto [S, runtime] = model.run(60, problem.T * problem.V * problem.L * problem.K, 23168);
    std::cout << S << '\n' << runtime << '\n';



    return 0;
}
