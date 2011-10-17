#pragma once

#include <iostream>
#include <vector>

namespace tsp {
namespace utils {

const double EPSILON = 1e-10;

std::vector<std::vector<double> > readData(std::istream &input) {
    int n;
    input >> n;
    std::vector<std::vector<double> > metric(n, std::vector<double>(n));
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            input >> metric[row][col];
        }
    }
    return metric;
}

}
}
