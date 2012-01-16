#include "branch_bound.h"
#include "utils.h"

#include <iostream>
#include <utility>
#include <vector>

#include <cstdlib>
#include <ctime>

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./bb <n>" << std::endl;
        return 1;
    }
    int n = atoi(argv[1]);
    srand(time(NULL));
    std::vector<std::pair<double, double> > points = tsp::utils::generateEuclidean(n);
    std::vector<std::vector<double> > metric = tsp::utils::getEuclideanMetric(points);
    tsp::core::BranchBound bb(metric);
    bb.run();
    return 0;
}