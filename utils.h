#pragma once

#include <iostream>
#include <utility>
#include <vector>

#include <cmath>
#include <cstdlib>

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

std::vector<std::pair<double, double> > generateEuclidean(int n) {
    std::vector<std::pair<double, double> > points;
    for (int i = 0; i < n; ++i) {
        double x = (rand() + 1.0) / (RAND_MAX + 0.0);
        double y = (rand() + 1.0) / (RAND_MAX + 0.0);
        points.push_back(std::make_pair(x, y));
    }
    return points;
}

template<typename T> T sqr(T x) {
    return x * x;
}

std::vector<std::vector<double> > getEuclideanMetric(const std::vector<std::pair<double, double> > &points) {
    int n = static_cast<int>(points.size());
    std::vector<std::vector<double> > metric(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            metric[i][j] = sqrt(sqr(points[i].first - points[j].first) + sqr(points[i].second - points[j].second));
        }
    }
    return metric;
}

}
}
