#include "utils.h"
#include "held_karp.h"
#include "heuristics.h"

#include <iostream>
#include <stdexcept>
#include <vector>

#include <cstddef>

int main() {
    try {
        std::cerr << "----------" << std::endl;
        std::vector<std::vector<double> > metric = tsp::utils::readData(std::cin);
        int n = static_cast<int>(metric.size());
        std::vector<std::vector<char> > emptyTour(n, std::vector<char>(n, -1));
        std::vector<std::vector<double> > fractionalTour(n, std::vector<double>(n));
        tsp::core::HeldKarpLowerBound heldKarp(metric, emptyTour);
        double lpBound = heldKarp.getBound(fractionalTour);
        std::cerr.precision(20);
        std::cerr << "HELD_KARP_BOUND " << lpBound << std::endl;
        int numFractionalVariables = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double curValue = fractionalTour[i][j];
                if (curValue > tsp::utils::EPSILON && curValue < 1.0 - tsp::utils::EPSILON) {
                    ++numFractionalVariables;
                }
            }
        }
        std::cerr << "FRACTIONAL_VARIABLES " << numFractionalVariables << std::endl;
        tsp::core::TSPHeuristics heuristicSolver(metric);
        std::vector<int> nearestNeighborTour = heuristicSolver.getNearestNeighborTour();
        std::cerr << "NEAREST_NEIGHBOR_TOUR " << tsp::utils::evaluateTour(metric, nearestNeighborTour) << std::endl;
        std::vector<int> greedyTour = heuristicSolver.getGreedyTour();
        std::cerr << "GREEDY_TOUR " << tsp::utils::evaluateTour(metric, greedyTour) << std::endl;
        std::vector<int> freshCopy(greedyTour);
        heuristicSolver.do2Opt(freshCopy);
        std::cerr << "GREEDY_2OPT_TOUR " << tsp::utils::evaluateTour(metric, freshCopy) << std::endl;
        freshCopy = greedyTour;
        heuristicSolver.doBest2Opt(freshCopy);
        std::cerr << "GREEDY_BEST2OPT_TOUR " << tsp::utils::evaluateTour(metric, freshCopy) << std::endl;
        freshCopy = greedyTour;
        heuristicSolver.do3Opt(freshCopy);
        std::cerr << "GREEDY_3OPT_TOUR " << tsp::utils::evaluateTour(metric, freshCopy) << std::endl;
        freshCopy = greedyTour;
        heuristicSolver.doBest3Opt(freshCopy);
        std::cerr << "GREEDY_BEST3OPT_TOUR " << tsp::utils::evaluateTour(metric, freshCopy) << std::endl;
    }
    catch (std::runtime_error &e) {
        std::cerr << "runtime error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
