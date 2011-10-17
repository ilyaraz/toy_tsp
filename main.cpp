#include "utils.h"
#include "held_karp.h"

#include <iostream>
#include <stdexcept>
#include <vector>

#include <cstddef>

int main() {
    try {
        std::vector<std::vector<double> > metric = tsp::utils::readData(std::cin);
        size_t n = metric.size();
        std::vector<std::vector<char> > emptyTour(n, std::vector<char>(n, -1));
        std::vector<std::vector<double> > fractionalTour(n, std::vector<double>(n));
        tsp::core::HeldKarpLowerBound heldKarp(metric, emptyTour);
        double lpBound = heldKarp.getBound(fractionalTour);
        std::cerr.precision(20);
        std::cerr << "Held-Karp bound: " << lpBound << std::endl;
    }
    catch (std::runtime_error &e) {
        std::cerr << "runtime error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
