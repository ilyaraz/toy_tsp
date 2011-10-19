#pragma once

#include "utils.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tsp {
namespace core {

class TSPHeuristics {
public:
    TSPHeuristics(const std::vector<std::vector<double> > &metric): metric(metric) {
    }

    std::vector<int> getNearestNeighborTour() {
        size_t n = metric.size();
        std::vector<int> tour;
        int startPoint = rand() % n;
        tour.push_back(startPoint);
        std::vector<char> was(n, 0);
        was[startPoint] = 1;
        while (tour.size() < n) {
            int lastPoint = tour[tour.size() - 1];
            int bestPoint = -1;
            for (size_t i = 0; i < n; ++i) {
                if (was[i]) {
                    continue;
                }
                if (bestPoint == -1 || metric[lastPoint][i] < metric[lastPoint][bestPoint]) {
                    bestPoint = i;
                }
            }
            tour.push_back(bestPoint);
            was[bestPoint] = 1;
        }
        return tour;
    }
    
    std::vector<int> getGreedyTour() {
        int n = static_cast<int>(metric.size());
        std::vector<std::pair<double, std::pair<int, int> > > edges(n * (n - 1) / 2);
        int numEdges = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                edges[numEdges] = std::make_pair(metric[i][j], std::make_pair(i, j));
                ++numEdges;
            }
        }
        std::sort(edges.begin(), edges.end());
        tsp::utils::DSU dsu(n);
        std::vector<int> degree(n, 0);
        int edgesAdded = 0;
        std::vector<int> neighbors(2 * n, -1);
        for (int i = 0; i < numEdges; ++i) {
            const std::pair<double, std::pair<int, int> > &edge = edges[i];
            int u = edge.second.first;
            int v = edge.second.second;
            if (degree[u] >= 2 || degree[v] >= 2) {
                continue;
            }
            if (edgesAdded != n - 1 && dsu.getRoot(u) == dsu.getRoot(v)) {
                continue;
            }
            dsu.connect(u, v);
            ++degree[u];
            ++degree[v];
            ++edgesAdded;
            addNeighbor(neighbors, u, v);
            addNeighbor(neighbors, v, u);
        }
        std::vector<char> was(n, 0);
        was[0] = 1;
        std::vector<int> tour(n);
        int tourLength = 1;
        tour[0] = 0;
        while (tourLength != n) {
            int lastPoint = tour[tourLength - 1];
            int nextPoint;
            if (!was[neighbors[2 * lastPoint]]) {
                nextPoint = neighbors[2 * lastPoint];
            }
            else {
                nextPoint = neighbors[2 * lastPoint + 1];
            }
            tour[tourLength++] = nextPoint;
            was[nextPoint] = 1;
        }
        return tour;
    }

    void do2Opt(std::vector<int> &tour) {
        int n = static_cast<int>(tour.size());
        int numRounds = 0;
        int numSwaps = 0;
        for (;;) {
            bool done = true;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    int c1 = tour[i];
                    int c2 = tour[(i + 1) % n];
                    int c3 = tour[j];
                    int c4 = tour[(j + 1) % n];
                    if (metric[c1][c2] + metric[c3][c4] > metric[c1][c3] + metric[c2][c4] + tsp::utils::EPSILON) {
                        ++numSwaps;
                        done = false;
                        int u = (i + 1) % n, v = j;
                        while (u < v) {
                            std::swap(tour[u], tour[v]);
                            ++u;
                            --v;
                        }
                    }
                }
            }
            if (done) {
                break;
            }
            ++numRounds;
        }
        std::cerr << "2OPT_ROUNDS " << numRounds << std::endl;
        std::cerr << "2OPT_SWAPS " << numSwaps << std::endl;
    }

    void doBest2Opt(std::vector<int> &tour) {
        int n = static_cast<int>(tour.size());
        int numSwaps = 0;
        for (;;) {
            int bestI = -1, bestJ = -1;
            double bestGain = 0.0;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    int c1 = tour[i];
                    int c2 = tour[(i + 1) % n];
                    int c3 = tour[j];
                    int c4 = tour[(j + 1) % n];
                    double gain = metric[c1][c2] + metric[c3][c4] - metric[c1][c3] - metric[c2][c4];  
                    if (gain > bestGain) {
                        bestGain = gain;
                        bestI = i;
                        bestJ = j;
                    }
                }
            }
            if (bestGain < tsp::utils::EPSILON) {
                break;
            }
            int u = (bestI + 1) % n, v = bestJ;
            while (u < v) {
                std::swap(tour[u], tour[v]);
                ++u;
                --v;
            }
            ++numSwaps;
        }
        std::cerr << "BEST2OPT_SWAPS " << numSwaps << std::endl;
    }
private:
    const std::vector<std::vector<double> > &metric;

    void addNeighbor(std::vector<int> &neighbors, int u, int v) {
        if (neighbors[2 * u] == -1) {
            neighbors[2 * u] = v;
            return;
        }
        if (neighbors[2 * u + 1] != -1) {
            throw std::runtime_error("already filled");
        }
        neighbors[2 * u + 1] = v;
    }
};

}
}
