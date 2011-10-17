#pragma once

#include <iostream>
#include <vector>

namespace tsp {
namespace utils {

std::vector<std::vector<double> > readData(std::istream &input) {
    size_t n;
    input >> n;
    std::vector<std::vector<double> > metric(n, std::vector<double>(n));
    for (size_t row = 0; row < n; ++row) {
        for (size_t col = 0; col < n; ++col) {
            input >> metric[row][col];
        }
    }
    return metric;
}

}
}
