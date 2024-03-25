#include <bits/stdc++.h>
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "genetic_algorithm.cpp"
#include "util/Turbine.cpp"
#include "util/Cable.cpp"
#include "SWFCR.cpp"

// 1: path, 2: filename.inf, 3: filename.turb, 4: filename.cbl
int main(int argc, char *argv[]) {
    assert(argc == 5);

    std::string dir = argv[1];

    read(dir + argv[2], dir + argv[3], dir + argv[4]);
    init_graph();


    Evolution ga(N);
    ga.run_evo(-1, 0.1, 10);



    return 0;
}
