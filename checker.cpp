#include <fstream>
#include <iostream>
#include <utility>
#include <set>
#include <vector>

#include "utils.h"

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "usage: ./checker <input-file> <output-file>" << std::endl;
        return 1; 
    }
    std::ifstream input(argv[1]);
    std::ifstream output(argv[2]);
    int n;
    input >> n;
    std::vector<std::pair<double, double> > points(n);
    for (int i = 0; i < n; ++i) {
        input >> points[i].first >> points[i].second;
    }
    std::vector<std::vector<double> > metric = tsp::utils::getEuclideanMetric(points);
    std::vector<int> answer(n);
    std::set<int> was;
    for (int i = 0; i < n; ++i) {
        output >> answer[i];
        if (answer[i] < 0 || answer[i] >= n) {
            std::cerr << "invalid node: " << answer[i] << std::endl;
            return 1;
        }
        if (was.count(answer[i])) {
            std::cerr << "not a permutation" << std::endl;
            return 1;
        }
        was.insert(answer[i]);
    }
    for (int i = 0; i < n; ++i) {
        std::cerr << answer[i] << " ";
    }
    std::cerr << std::endl;
    double weight = 0.0;
    for (int i = 0; i < n; ++i) {
        weight += metric[answer[i]][answer[(i + 1) % answer.size()]];
    }
    std::cerr << "weight: " << weight << std::endl;
    return 0;
}
