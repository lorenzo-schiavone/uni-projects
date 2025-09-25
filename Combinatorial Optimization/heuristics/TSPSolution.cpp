#include "TSPSolution.h"
#include <numeric>
#include <random>


TSPSolution::TSPSolution() : cost(0) {}

TSPSolution::TSPSolution(const TSPInstance& instance) : cost(0) {
    tour.resize(instance.numCities);
    std::iota(tour.begin(), tour.end(), 0); // 0123456789...
    tour.push_back(0);
}

TSPSolution& TSPSolution::operator=(const TSPSolution& right) {
    if (this == &right) return *this;
    tour = right.tour;
    cost = right.cost;
    return *this;
}

double evaluateTour(const TSPInstance& instance, const TSPSolution& sol) {
    double total = 0.0;
    for (int i = 0; i < instance.numCities; i++) {
        total += instance.costMatrix[sol.tour[i]][sol.tour[i+1]];
    }
    return total;
}

// CREATE SOLUTION:
// RANDOM
TSPSolution createRandomSolution(const TSPInstance& instance, std::mt19937 &rng) {
    TSPSolution sol(instance);
    std::shuffle(sol.tour.begin() + 1, sol.tour.end()-1, rng);
    sol.cost = evaluateTour(instance, sol);
    return sol;
}

// HEURISTIC
TSPSolution farthestInsertion (const TSPInstance& tsp){
    int numCity = tsp.numCities;
    TSPSolution sol(tsp);
    std::vector<bool> unvisited(numCity, true);
    unvisited[0]=false;

    int i_best = -1;
    double max_dist = 0;
    for ( int i = 1 ; i < numCity-1; ++i ) {
            if (tsp.costMatrix[0][i] > max_dist) {
                max_dist = tsp.costMatrix[0][i];
                i_best = i;
            }
    }
    unvisited[i_best] = false;
    std::vector<int> tour = {0, i_best, 0};
    while (tour.size()<numCity+1) {
        int r = -1;
        double best_value = 0;

        for (int v=1; v<numCity; v++) {
            if (!unvisited[v])continue;
            double max_dist_to_tour = 0;
            for (int j=0; j<tour.size()-1; j++) {
                int w=tour[j];
                max_dist_to_tour = std::max(max_dist_to_tour, tsp.costMatrix[v][w]);
            }
            if (max_dist_to_tour > best_value) {
                best_value = max_dist_to_tour;
                r = v;
            }
        }

        int best_pos = -1;
        double best_cost = tsp.infinite;

        for (int idx = 0; idx < tour.size() - 1; idx++) {
            int i = tour[idx], j = tour[idx + 1];
            double insertion_cost = tsp.costMatrix[i][r] + tsp.costMatrix[r][j] - tsp.costMatrix[i][j];

            if (insertion_cost < best_cost) {
                best_cost = insertion_cost;
                best_pos = idx + 1;
            }
        }
        tour.insert(tour.begin() + best_pos, r);
        unvisited[r]=false;
    }
    sol.tour = tour; 
    sol.cost = evaluateTour(tsp,sol);
    return sol;
}

// nearestInsertion
TSPSolution nearestInsertion (const TSPInstance& tsp) {
    int numCity = tsp.numCities;
    TSPSolution sol(tsp);

    std::vector<bool> unvisited(numCity, true);
    unvisited[0] = false;

    int i_best = -1;
    double min_dist = tsp.infinite;
    for (int i = 1; i < numCity; ++i) {
        if (tsp.costMatrix[0][i] < min_dist) {
            min_dist = tsp.costMatrix[0][i];
            i_best = i;
        }
    }

    unvisited[i_best] = false;
    std::vector<int> tour = {0, i_best, 0};

    while (static_cast<int>(tour.size()) < numCity + 1) {
        int r = -1;
        double best_value = tsp.infinite;

        for (int v = 1; v < numCity; ++v) {
            if (!unvisited[v]) continue;
            double min_dist_to_tour = tsp.infinite;
            for (int j = 0; j < static_cast<int>(tour.size()) - 1; ++j) {
                int w = tour[j];
                min_dist_to_tour = std::min(min_dist_to_tour, tsp.costMatrix[v][w]);
            }
            if (min_dist_to_tour < best_value) {
                best_value = min_dist_to_tour;
                r = v;
            }
        }

        int best_pos = -1;
        double best_cost = tsp.infinite;
        for (int idx = 0; idx < static_cast<int>(tour.size()) - 1; ++idx) {
            int i = tour[idx], j = tour[idx + 1];
            double insertion_cost = tsp.costMatrix[i][r] + tsp.costMatrix[r][j] - tsp.costMatrix[i][j];
            if (insertion_cost < best_cost) {
                best_cost = insertion_cost;
                best_pos = idx + 1;
            }
        }
        tour.insert(tour.begin() + best_pos, r);
        unvisited[r] = false;
    }

    sol.tour = tour;
    sol.cost = evaluateTour(tsp, sol);
    return sol;
}
