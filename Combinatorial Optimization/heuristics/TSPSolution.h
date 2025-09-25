#ifndef TSP_SOLUTION_HPP
#define TSP_SOLUTION_HPP

#include "TSP.h"
#include <vector>
#include <random>

struct TSPSolution {
    std::vector<int> tour;
    double cost;
    TSPSolution();
    explicit TSPSolution(const TSPInstance& instance);
    TSPSolution& operator=(const TSPSolution& right);
};

double evaluateTour(const TSPInstance& instance, const TSPSolution& sol);
TSPSolution createRandomSolution(const TSPInstance& instance, std::mt19937 &rng);
TSPSolution farthestInsertion (const TSPInstance& tsp);
TSPSolution nearestInsertion (const TSPInstance& tsp);

#endif // TSP_SOLUTION_HPP