#include "utilities.h"
#include "genetic_algorithm.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>

int main (int argc, char const *argv[]){
    if (argc < 2) { 
        std::cerr << "Usage: " << argv[0] << " <instance>" << std::endl;
        return 1;
    }
    bool verbose = false;
    if (argc > 2) {
        if (strcmp(argv[2], "-v") == 0) {
            verbose = true;
        }
    }

    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::chrono::high_resolution_clock::time_point startTime;
    
    std::string base_path1 = "../data/cost_matrix/cost";
    std::string extension1 = ".dat";
    std::string path1 = base_path1 + argv[1] + extension1;
    const char* costs_path = path1.c_str();

    std::string base_path2 = "../solutions/localSearch/sol";
    std::string extension2 = ".txt";
    std::string path2 = base_path2 + argv[1] + extension2;
    const char* sol_path = path2.c_str();

    std::string base_path3 = "../results/localSearch_summary.csv";
    const char* save_path = base_path3.c_str();

    TSPInstance tspInstance;
    tspInstance.read(costs_path);
    
    double TimeLimit = getTimeLimit("./parameters/timeLimit.txt", tspInstance.numCities);
    
    int LScount = 0;
        
    startTime = std::chrono::high_resolution_clock::now();
    TSPSolution bestSol = farthestInsertion(tspInstance);
        
    localSearch(tspInstance, bestSol, verbose);
    double currTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();
    double bestVal = bestSol.cost;
    std::cout << "Objval: " << bestVal << ", " << "in " << currTime << " sec\n"; // 
    double lastImprovTime=currTime;

    while (currTime<TimeLimit){
        if (LScount > 2){
            verbose = false;
        }
        if (verbose){
            std::cout << "\nrestart LS";
        }
        TSPSolution sol = createRandomSolution(tspInstance, rng);
        localSearch(tspInstance, sol, verbose);
        currTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();

        if (sol.cost<bestVal-1e-14){
            bestSol = sol;
            bestVal = sol.cost;
            currTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();
            std::cout << "New Objval: " << bestVal << ", " << "in " << currTime << " sec\n"; // 
            lastImprovTime = currTime;

        }
        currTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();
        LScount++;
    }

    writeSol(sol_path, bestSol.tour, bestVal);
    writeTimeObj(save_path, argv[1], (double) bestSol.cost, (double) lastImprovTime);
    if (verbose){
        std::cout << "LS done: " << LScount << std::endl;
    }
    std::cout << "Best Objective value: " << bestVal  << "\n";

    return 0;
}