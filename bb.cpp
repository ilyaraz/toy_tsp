#include "branch_bound.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <utility>
#include <string>
#include <vector>

#include <cstdlib>
#include <ctime>

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: ./bb -n <n>" << std::endl;
        std::cerr << "       ./bb -f <file>" << std::endl;
        return 1;
    }
    std::string mode(argv[1]);
    if (mode == "-n") {
        int n = atoi(argv[2]);
        srand(time(NULL));
        std::vector<std::pair<double, double> > points = tsp::utils::generateEuclidean(n);
        std::vector<std::vector<double> > metric = tsp::utils::getEuclideanMetric(points);
        tsp::core::BranchBound bb(metric);
        bb.run();
        std::ofstream output("points.txt");
        output << n << std::endl;
        output.setf(std::ios::fixed);
        output.precision(20);
        for (int i = 0; i < n; ++i) {
            output << points[i].first << " " << points[i].second << std::endl;
        }
    }
    else if (mode == "-f") {
        std::ifstream input(argv[2]);
        int n;
        input >> n;
        std::vector<std::pair<double, double> > points(n);
        for (int i = 0; i < n; ++i) {
            input >> points[i].first >> points[i].second;
        }
        std::vector<std::vector<double> > metric = tsp::utils::getEuclideanMetric(points);
        tsp::core::BranchBound bb(metric);
        bb.run();
    }
    else {
        std::cerr << "Usage: ./bb -n <n>" << std::endl;
        std::cerr << "       ./bb -f <file>" << std::endl;
        return 1;
    }
    return 0;
}
