#include "utils.h"

#include <iostream>
#include <vector>

#include <cstddef>

int main() {
    std::vector<std::vector<double> > metric = tsp::utils::readData(std::cin);
    size_t n = metric.size();
    std::cerr.precision(20);
    for (size_t row = 0; row < n; ++row) {
        for (size_t col = 0; col < n; ++col) {
            std::cerr << metric[row][col] << " ";
        }
        std::cerr << std::endl;
    }
    return 0;
}
