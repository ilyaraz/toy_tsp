#pragma once

#include <iostream>
#include <utility>
#include <vector>

#include <cmath>
#include <cstdlib>

namespace tsp {
namespace utils {

const double EPSILON = 1e-10;
const double PI = acos(-1.0);

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

double getNormal() {
    double u = (rand() + 1.0) / (RAND_MAX + 0.0);
    double v = (rand() + 1.0) / (RAND_MAX + 0.0);
    return sqrt(-2.0 * log(u)) * cos(2 * PI * v);
}

std::vector<std::pair<double, double> > generateNormal(int n) {
    std::vector<std::pair<double, double> > points;
    for (int i = 0; i < n; ++i) {
        double x = 0.5 + 0.1 * getNormal(); 
        double y = 0.5 + 0.1 * getNormal();
        points.push_back(std::make_pair(x, y));
    }
    return points;
}

template<typename T> T sqr(T x) {
    return x * x;
}

template<typename T> T abs(T x) {
    return (x > 0) ? x : -x;
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

double evaluateTour(const std::vector<std::vector<double> > &metric, const std::vector<int> &tour) {
    double length = 0.0;
    int n = static_cast<int>(tour.size());
    for (int i = 0; i < n; ++i) {
        length += metric[tour[i]][tour[(i + 1) % n]];
    }
    return length;
}

class DSU {
public:
    DSU(int n): parent(n) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    int getRoot(int x) {
        if (x != parent[x]) {
            parent[x] = getRoot(parent[x]);
        }
        return parent[x];
    }

    void connect(int x, int y) {
        link(getRoot(x), getRoot(y));
    }
private:
    std::vector<int> parent;

    void link(int x, int y) {
        if (rand() % 2) {
            parent[x] = y;
        }
        else {
            parent[y] = x;
        }
    }
};

}
}
