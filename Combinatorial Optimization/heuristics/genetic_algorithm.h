#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "TSP.h"
#include "TSPSolution.h"
#include "localSearch.h"
#include "utilities.h"
#include <vector>
#include <random>


TSPSolution runGeneticAlgorithm(TSPInstance tspInstance, int popSize, int newPopSize, double lsRate, double MAXmutationRate, int MAXIT, double replace_factor, std::mt19937& rng, bool verbose=false);
#endif // INITIALIZATION_H