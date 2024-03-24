#include <bits/stdc++.h>
#include "util/DebugTemplate.cpp"
#include "util/Constants.cpp"
#include "genetic_algorithm.cpp"
#include "util/Turbine.cpp"
#include "util/Cable.cpp"
#include "SWFCR.cpp"

// 1: path, 2: filename.turb, 3: filename.cbl
int main(int argc, char *argv[]) {
    assert(argc == 4);

    std::string dir = argv[1];

    read(dir + "\\" + argv[2], dir + "\\" + argv[3]);
    init_graph();


    Evolution ga;


    return 0;
}
