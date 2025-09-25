#include "TSP.h"
#include "TSPSolution.h"
#include "localSearch.h"
#include "utilities.h"
#include <random>
#include <iostream>
#include <cmath>

// Global Variables
std::uniform_real_distribution<> prob01 (0.0,1.0);
std::uniform_int_distribution<> randomIndex;
int numCity;

// ----------------------------------------------------------------------
// POPULATION OPERATION

std::vector<TSPSolution> initPopulation(const TSPInstance& instance, int popSize, std::mt19937 &rng) {
    std::vector<TSPSolution> population;
    population.reserve(popSize);
    population.push_back(farthestInsertion(instance));
    population.push_back(nearestInsertion(instance));

    TSPSolution sol= createRandomSolution(instance, rng);
    localSearch2opt(instance, sol);
    population.push_back(sol);

    int leftPop = popSize - population.size();
    for (int i=0; i<leftPop; i++ ) {
        population.push_back(createRandomSolution(instance, rng));
    }
    return population;
}

// probability for linear ranking
std::vector<double> linearRanking(const int N){ 
    std::vector<double> rankProbabilities(N);
    double cumsum=0.0;
    for (int i=0; i<N; i++){
        cumsum = cumsum + (2.0 *(i+1)) /(N*(N+1));
        rankProbabilities[i] = cumsum;
    }
    return rankProbabilities;
}

// linear selection from population of selectedSize using probabilities from rankProbabilities !! here it sorts and returns unsor
std::vector<TSPSolution> linearSelection(const std::vector<TSPSolution>& population, const int selectedSize, std::vector<double>& rankProbabilities, std::mt19937 &rng) {
    int size = population.size();

    std::vector<int> indices(size);
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
        [&](int a, int b) {
            return population[a].cost > population[b].cost;
        });

    std::vector<bool> unselectedIndeces(size, true); // avoid multiple selection
    std::vector<TSPSolution> newPopulation(selectedSize);
    int cc = 0;
    while (cc < selectedSize){
        double p = prob01(rng);
        for (int i = 0; i < size; ++i) {
            if (p <= rankProbabilities[i] || i == size-1) {
                int idx = indices[i];
                if (unselectedIndeces[idx]){
                    unselectedIndeces[idx] = false;
                    newPopulation[cc] = population[idx];
                    cc++;
                }
                break;
            }
        }
    }
    return newPopulation;
}

// CROSSOVER
TSPSolution crossover(const TSPSolution& parent1, const TSPSolution& parent2, std::mt19937 &rng) {
    
    TSPSolution child = parent1;
    int start = randomIndex(rng);
    int end = randomIndex(rng);
    if (start > end) std::swap(start, end);
    
    std::vector<bool> unvisited(numCity, false);
    for (int i=start; i<= end; i++){
        unvisited[parent1.tour[i]]=true;
    }

    int current=start;
    for (int i=1; i<numCity;i++){
        if (unvisited[parent2.tour[i]]){
            child.tour[current] = parent2.tour[i];
            unvisited[parent2.tour[i]]=false;
            current++;
        }
    }

    return child;
}

// NATURAL SELECTION
std::vector<TSPSolution> naturalSelection( std::vector<TSPSolution>& population, int popSize, std::vector<double>& rankProbabilities, std::mt19937 &rng){
    const int elite = popSize / 20; // 5%

    std::partial_sort(population.begin(), population.begin() + elite, population.end(),
        [](const TSPSolution& a, const TSPSolution& b){
            return a.cost < b.cost;
        }); // find elite

    std::vector<TSPSolution> selected;
    selected.reserve(popSize);
    selected.insert(selected.end(), population.begin(), population.begin() + elite); // insert elite
    std::vector<TSPSolution> remaining = linearSelection({ population.begin() + elite, population.end() }, popSize - elite, rankProbabilities, rng);
    
    selected.insert(selected.end(), remaining.begin(), remaining.end());

    return selected;
}

// MUTATION: swap
void mutate(TSPSolution &sol, double mutationRate, std::mt19937 &rng) {
    if (prob01(rng) < mutationRate) {
        int i = randomIndex(rng);
        int j = randomIndex(rng);
        std::swap(sol.tour[i], sol.tour[j]);
    }
}

std::vector<double> temperatureSchedule(int N, double MAXmutationRATIO) {
    std::vector<double> mutationRatio(N);
    for(int iter=0; iter<N; iter++) {
        double progress = static_cast<double>(iter)/N;
        mutationRatio[iter] = (MAXmutationRATIO * (1.0 - progress));  
    }
    return mutationRatio;
}

// METRIC ON DIVERSITY
double hammingDistance(const TSPSolution& sol1,const TSPSolution& sol2 ){
    int diff1 = 0;
    int diff2 = 0;
    for (int k = 1; k < numCity; k++) {
        if (sol1.tour[k] != sol2.tour[k]) // 1 if different
            diff1++;
        if (sol1.tour[k] != sol2.tour[numCity-k])
            diff2++;
    }
    return std::min(diff1,diff2)/static_cast<double>(numCity-1);
}

double averageHammingDistance(const std::vector<TSPSolution>& population) {
    int n = population.size();
    double total = 0.0;
    int count = 0;
    for (int i = 0; i < n-1; i++) {
        for (int j = i+1; j < n; j++) {
            total += hammingDistance(population[i],population[j]);
            count++;
        }
    }
    return total / count;
}

// ----------------------------------------------------------------------
// ALGORITHM

TSPSolution runGeneticAlgorithm(TSPInstance tspInstance, int popSize, int newPopSize,  double lsRate, double MAXmutationRate, int MAXIT, double replace_factor, std::mt19937& rng, bool verbose){
    if (verbose){
            std::cout << "popSize: " << popSize 
            << ", newPopSize: " << newPopSize 
            << ", maxit: " << MAXIT << std::endl;
    }

    std::uniform_int_distribution<> rndpop (1, popSize-1); // non prende il primo 
    numCity = tspInstance.numCities;
    randomIndex.param(std::uniform_int_distribution<>::param_type(1, numCity-1)); // load correct size for choice random index in a tour: for mutate and ordercrossover
    std::vector<double> rankPopulation = linearRanking(popSize);
    double avgHamming;
    int matingPoolSize = std::ceil((1+std::sqrt(1+8*newPopSize))/2); // smallest n st (matingPoolSize+1)*matingPoolSize / 2 > newPopSize
    int effectiveNewPopSize = (matingPoolSize+1)*(matingPoolSize)/2;
    std::vector<double> rankNaturalSelection = linearRanking(popSize*.95 + effectiveNewPopSize);

    std::vector<TSPSolution> population = initPopulation(tspInstance, popSize, rng);
    TSPSolution bestSol = *min_element(population.begin(), population.end(), [](const TSPSolution& a, const TSPSolution& b) {
        return a.cost < b.cost;
    });

    std::vector<double> mutationRatio = temperatureSchedule(MAXIT, MAXmutationRate);

    int resetTemp=0; int resetCounter = 0;
    int replacePop = std::floor(popSize*replace_factor);
    for (int iter=0; iter<MAXIT; iter++) { 

        std::vector<TSPSolution> matingPool = linearSelection(population, matingPoolSize, rankPopulation, rng);
                
        /// MATING ALL THE MATING POOL
        for (int i=0;i<matingPool.size()-1;i++){
            for (int j=i+1;j<matingPool.size();j++){
                // CHILD GENERATION
                TSPSolution child = crossover(matingPool[i],matingPool[j], rng);

                mutate(child, mutationRatio[iter-resetTemp], rng);

                if (prob01(rng) < lsRate){ // local search on very fex child
                    localSearch2opt(tspInstance, child);
                }
                child.cost = evaluateTour(tspInstance, child);
                if (child.cost < bestSol.cost-1e-14){
                    bestSol = child; 
                }
                population.push_back(child);
            }
        }
        population = naturalSelection(population, popSize, rankNaturalSelection, rng); // elite + linear ranking

        avgHamming = averageHammingDistance(population);
        if (avgHamming<.15){ //population management - injection of new random
            if (replacePop < popSize/2){
                std::vector<bool> notReplacedYet(popSize, true);
                for (int i= 0; i<replacePop;i++){
                    int idx;
                    do {
                        idx = rndpop(rng);
                    }while (!notReplacedYet[idx]);
                    population[idx] = createRandomSolution(tspInstance, rng);
                    notReplacedYet[idx] = false; 
                }
            }
            else {
                std::vector<int> idx(population.size()-1);
                std::iota(idx.begin(), idx.end(), 1); // preserve best
                std::shuffle(idx.begin(), idx.end(), rng);
                for (int k = 0; k < replacePop; k++) {
                    population[idx[k]] = createRandomSolution(tspInstance, rng);
                }
            }
            resetTemp = iter; // bump again temperature
            resetCounter++;
            mutationRatio = temperatureSchedule(MAXIT - resetTemp, MAXmutationRate);
        }

        if (verbose){
            if (iter%100 == 0){
                std::cout << "Iteration " << iter 
                << " | Best Tour Length: " << bestSol.cost 
                << " | Avg Hamming: " << avgHamming << std::endl;
                }
        }

    }
    // last minute improvement
    localSearch(tspInstance, bestSol);
    if (verbose){
        std::cout << "Population injection count: " << resetCounter << "\n";
    }
    return bestSol;
}