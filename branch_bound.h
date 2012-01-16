#pragma once

#include "held_karp.h"
#include "heuristics.h"
#include "utils.h"

#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include <cstdlib>

namespace tsp { 
namespace core {

struct BranchBoundGlobalState {
    double upperBound;
    std::vector<int> bestTour;

    BranchBoundGlobalState(double upperBound, const std::vector<int> &bestTour):
        upperBound(upperBound),
        bestTour(bestTour) {}
};

struct BranchBoundNode {
    double lowerBound;
    std::vector<std::vector<char> > partialTour;
    std::vector<std::vector<double> > bestFractionalTour;

    BranchBoundNode(double lowerBound,
                    const std::vector<std::vector<char> > &partialTour,
                    const std::vector<std::vector<double> > &bestFractionalTour):
        lowerBound(lowerBound),
        partialTour(partialTour),
        bestFractionalTour(bestFractionalTour) {}

    bool isIntegral() const {
        int n(bestFractionalTour.size());
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double tourIJ = bestFractionalTour[i][j];
                if (tourIJ > tsp::utils::EPSILON && tourIJ < 1.0 - tsp::utils::EPSILON) {
                    return false;
                }
            }
        }
        return true;
    }
};

struct BranchBoundNodesComparator {
    bool operator()(const BranchBoundNode &node1, const BranchBoundNode &node2) const {
        return node1.lowerBound < node2.lowerBound;
    }
};

class BranchBound {
    public:
        BranchBound(const std::vector<std::vector<double> > &metric): n(metric.size()), metric(metric) {
        }
        
        void run() {
            TSPHeuristics solver(metric);
            std::vector<int> heuristicTour = solver.getGreedyTour();
            solver.doBest3Opt(heuristicTour);
            BranchBoundGlobalState globalState = createBranchBoundGlobalState(heuristicTour);
            std::vector<std::vector<char> > emptyPartialTour(n, std::vector<char>(n, -1));
            std::set<BranchBoundNode, BranchBoundNodesComparator> branchBoundNodes;
            branchBoundNodes.insert(createBranchBoundNode(emptyPartialTour));
            for (;;) {
                for (std::set<BranchBoundNode, BranchBoundNodesComparator>::const_iterator it = branchBoundNodes.begin();
                        it != branchBoundNodes.end();
                        ++it) {
                    std::cerr << it->lowerBound << " ";
                }
                std::cerr << std::endl;
                if (branchBoundNodes.empty()) {
                    std::cerr << "Pruning..." << std::endl;
                    std::cerr << "Answer: " << globalState.upperBound << std::endl;
                    break;
                }
                BranchBoundNode currentNode = *branchBoundNodes.begin();
                branchBoundNodes.erase(branchBoundNodes.begin());
                std::cerr << currentNode.lowerBound << " vs " << globalState.upperBound << std::endl;
                std::vector<BranchBoundNode> children = getChildren(currentNode, globalState);
                branchBoundNodes.insert(children.begin(), children.end());
            }
        }
    private:
        int n;
        const std::vector<std::vector<double> > &metric;

        BranchBoundGlobalState createBranchBoundGlobalState(const std::vector<int> &tour) {
            double cost = tsp::utils::evaluateTour(metric, tour);
            return BranchBoundGlobalState(cost, tour);
        }

        BranchBoundNode createBranchBoundNode(const std::vector<std::vector<char> > &partialTour) {
            HeldKarpLowerBound solver(metric, partialTour);
            std::vector<std::vector<double> > bestFractionalTour(n, std::vector<double>(n));
            double lowerBound = solver.getBound(bestFractionalTour);
            return BranchBoundNode(lowerBound, partialTour, bestFractionalTour);
        }

        std::vector<BranchBoundNode> getChildren(const BranchBoundNode &currentNode,
                                                 BranchBoundGlobalState &globalState) {
            if (currentNode.lowerBound > globalState.upperBound - tsp::utils::EPSILON) {
                return std::vector<BranchBoundNode>();
            }
            if (currentNode.isIntegral()) {
                std::cerr << "update: " << globalState.upperBound << " -> " << currentNode.lowerBound << std::endl;
                int tmp;
                std::cin >> tmp;
                globalState.upperBound = currentNode.lowerBound;
                return std::vector<BranchBoundNode>();
            }
            /*Applying a heuristic*/
            TSPHeuristics solver(metric);
            std::vector<int> heuristicTour = solver.getGreedyTour(currentNode.bestFractionalTour);
            solver.doBest3Opt(heuristicTour);
            double cost = tsp::utils::evaluateTour(metric, heuristicTour);
            if (cost < globalState.upperBound - tsp::utils::EPSILON) {
                std::cerr << "heuristic update: " << globalState.upperBound << " -> " << cost << std::endl;
                int tmp;
                std::cin >> tmp;
                globalState.upperBound = cost;
            }
            std::vector<std::vector<double> > _partialTour(n, std::vector<double>(n));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    _partialTour[i][j] = currentNode.partialTour[i][j];
                }
            }
            heuristicTour = solver.getGreedyTour(_partialTour);
            solver.doBest3Opt(heuristicTour);
            cost = tsp::utils::evaluateTour(metric, heuristicTour);
            if (cost < globalState.upperBound - tsp::utils::EPSILON) {
                std::cerr << "heuristic update 2: " << globalState.upperBound << " -> " << cost << std::endl;
                int tmp;
                std::cin >> tmp;
                globalState.upperBound = cost;
            }
            /**/
            double bestDifference = 10.0;
            std::pair<int, int> branchPair(-1, -1);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    double difference = tsp::utils::abs(0.5 - currentNode.bestFractionalTour[i][j]);
                    if (difference < bestDifference) {
                        bestDifference = difference;
                        branchPair = std::make_pair(i, j);
                    }
                }
            }
            std::cerr << "branching on " << currentNode.bestFractionalTour[branchPair.first][branchPair.second] << std::endl;
            std::vector<BranchBoundNode> children;
            for (int i = 0; i < 2; ++i) {
                std::vector<std::vector<char> > partialTour(currentNode.partialTour);
                partialTour[branchPair.first][branchPair.second] = partialTour[branchPair.second][branchPair.first] = i;
                BranchBoundNode child = createBranchBoundNode(partialTour);
                if (child.lowerBound > globalState.upperBound - tsp::utils::EPSILON) {
                    continue;
                }
                if (child.isIntegral()) {
                    std::cerr << "update child: " << globalState.upperBound << " -> " << child.lowerBound << std::endl;
                    int tmp;
                    std::cin >> tmp;
                    globalState.upperBound = child.lowerBound;
                    continue;
                }
                children.push_back(child);
            }
            return children;
        }
};

}
}
