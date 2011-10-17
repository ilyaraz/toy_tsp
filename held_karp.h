#pragma once

#include "utils.h"

extern "C" {
#include "gurobi_c.h"
}

#include <lemon/list_graph.h>
#include <lemon/hao_orlin.h>

#include <stdexcept>
#include <vector>

namespace tsp {
namespace core {

class HeldKarpLowerBound {
    public:
        HeldKarpLowerBound(const std::vector<std::vector<double> > &metric,
                           const std::vector<std::vector<char> > &partialTour):
            metric(metric),
            partialTour(partialTour) {
        }

        double getBound(std::vector<std::vector<double> > &fractionalTour) {
            int n = static_cast<int>(metric.size());
            GRBenv *gurobiEnvironment;
            int result = GRBloadenv(&gurobiEnvironment, NULL);
            if (result != 0) {
                throw std::runtime_error("Can't create Gurobi environment");
            }
            result = GRBsetintparam(gurobiEnvironment, "OutputFlag", 0);
            if (result != 0) {
                throw std::runtime_error("Can't turn off Gurobi output");
            }
            int numVars = n * (n - 1) / 2;
            int counter = 0;
            std::vector<double> objectiveCoefficient(numVars);
            std::vector<double> lowerBound(numVars, 0.0);
            std::vector<double> upperBound(numVars, 1.0);
            std::vector<std::vector<int> > variableID(n, std::vector<int>(n, -1));
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    variableID[i][j] = variableID[j][i] = counter;
                    objectiveCoefficient[counter] = metric[i][j];
                    ++counter;
                }
            }
            GRBmodel *gurobiModel;
            result = GRBnewmodel(gurobiEnvironment,
                                 &gurobiModel,
                                 "HeldKarpLP",
                                 numVars,
                                 &objectiveCoefficient[0],
                                 &lowerBound[0],
                                 &upperBound[0],
                                 NULL,
                                 NULL);
            if (result != 0) {
                throw std::runtime_error("Can't create Gurobi model");
            }
            int numNonZeroes = n * (n - 1);
            std::vector<int> cbeg(n);
            std::vector<int> cind(numNonZeroes);
            std::vector<double> cval(numNonZeroes, 1.0);
            std::vector<char> sense(n, GRB_EQUAL);
            std::vector<double> rhs(n, 2.0);
            for (int i = 0; i < n; ++i) {
                cbeg[i] = i * (n - 1);
            }
            for (int i = 0; i < n; ++i) {
                counter = 0;
                for (int j = 0; j < n; ++j) {
                    if (i == j) {
                        continue;
                    }
                    int entryID = i * (n - 1) + counter;
                    cind[entryID] = variableID[i][j];
                    ++counter;
                }
            }
            result = GRBaddconstrs(gurobiModel, n, numNonZeroes, &cbeg[0], &cind[0], &cval[0], &sense[0], &rhs[0], NULL);
            if (result != 0) {
                throw std::runtime_error("Can't add constraints to Gurobi model");
            }
            result = GRBupdatemodel(gurobiModel);
            if (result != 0) {
                throw std::runtime_error("Can't update Gurobi model");
            }
            int numCuttingPlanes = 0;
            for (;;) {
                result = GRBoptimize(gurobiModel);
                if (result != 0) {
                    throw std::runtime_error("Can't optimize Gurobi model");
                }
                if (!numCuttingPlanes) {
                    double bimatchingBound;
                    result = GRBgetdblattr(gurobiModel, "ObjVal", &bimatchingBound);
                    if (result != 0) {
                        throw std::runtime_error("Can't query objective function");
                    }
                    std::cerr << "2-matching bound: " << bimatchingBound << std::endl;
                }
                std::vector<double> vars(numVars);
                result = GRBgetdblattrarray(gurobiModel, "X", 0, numVars, &vars[0]);
                std::vector<char> cuttingPlane(n);
                bool invalidSolution = getCuttingPlane(n, vars, cuttingPlane);
                if (!invalidSolution) {
                    break;
                }
                ++numCuttingPlanes;
                numNonZeroes = 0;
                for (int i = 0; i < n; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        if (cuttingPlane[i] != cuttingPlane[j]) {
                            cind[numNonZeroes] = variableID[i][j];
                            ++numNonZeroes;
                        }
                    }
                }
                result = GRBaddconstr(gurobiModel, numNonZeroes, &cind[0], &cval[0], GRB_GREATER_EQUAL, 2.0, NULL);
                if (result != 0) {
                    throw std::runtime_error("Can't add a constraint to Gurobi model");
                }
                result = GRBupdatemodel(gurobiModel);
                if (result != 0) {
                    throw std::runtime_error("Can't update Gurobi model");
                }
            }
            std::cerr << numCuttingPlanes << " cutting planes" << std::endl;
            double hkValue;
            result = GRBgetdblattr(gurobiModel, "ObjVal", &hkValue);
            if (result != 0) {
                throw std::runtime_error("Can't query objective function");
            }
            result = GRBfreemodel(gurobiModel);
            if (result != 0) {
                throw std::runtime_error("Can't free Gurobi model");
            }
            GRBfreeenv(gurobiEnvironment);
            return hkValue;
        }
    private:
        const std::vector<std::vector<double> > &metric;
        const std::vector<std::vector<char> > &partialTour;

        bool getCuttingPlane(int n, const std::vector<double> &vars, std::vector<char> &cuttingPlane) {
            int counter = 0;
            lemon::ListGraph network;
            std::vector<lemon::ListGraph::Node> nodes(n);
            for (int i = 0; i < n; ++i) {
                nodes[i] = network.addNode();
            }
            lemon::ListGraph::EdgeMap<double> capacities(network);
            int numEdges = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    if (vars[counter] > tsp::utils::EPSILON) {
                        ++numEdges;
                        lemon::ListGraph::Edge edge = network.addEdge(nodes[i], nodes[j]);
                        capacities[edge] = vars[counter];
                    }
                    ++counter;
                }
            }
            lemon::HaoOrlin<lemon::ListGraph, lemon::ListGraph::EdgeMap<double> > cutSolver(network, capacities);
            cutSolver.init(nodes[0]);
            cutSolver.calculateOut();
            if (cutSolver.minCutValue() > 2.0 - tsp::utils::EPSILON) {
                return false;
            }
            lemon::ListGraph::NodeMap<bool> minCut(network);
            cutSolver.minCutMap(minCut);
            for (int i = 0; i < n; ++i) {
                cuttingPlane[i] = minCut[nodes[i]];
            }
            return true;
        }
};

}
}
