#include "utilities.h"
#include "genetic_algorithm.h"
#include "TSP.h"
#include <iostream>
#include <chrono>
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

    // paths for read cost matrix and write solution and result
    std::string base_path1 = "../data/cost_matrix/cost";
    std::string extension1 = ".dat";
    std::string path1 = base_path1 + argv[1] + extension1;
    const char* costs_path = path1.c_str();

    std::string base_path2 = "../solutions/genetic/sol";
    std::string extension2 = ".txt";
    std::string path2 = base_path2 + argv[1] + extension2;
    const char* sol_path = path2.c_str();

    std::string base_path3 = "../results/genetic_summary.csv";
    const char* save_path = base_path3.c_str();

    std::chrono::high_resolution_clock::time_point startTime;
    startTime = std::chrono::high_resolution_clock::now();
    
    TSPInstance tspInstance;
    tspInstance.read(costs_path);
    // get parameters
    std::string param_file = "./parameters/genetic.txt";
    std::tuple<int, int, double, double, double, int> param = getGeneticParam(param_file.c_str(), tspInstance.numCities);
    int popSize = std::get<0>(param);       
    int newPopSize = std::get<1>(param);    
    double lsRate =std::get<2>(param);      
    double MAXmutationRate = std::get<3>(param);
    double replace_factor = std::get<4>(param);
    int MAXIT = std::get<5>(param);
    
    TSPSolution bestSol = runGeneticAlgorithm(tspInstance, popSize, newPopSize, lsRate, MAXmutationRate, MAXIT, replace_factor, rng, verbose);
    double currTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();

    std::cout << "Objval: " << bestSol.cost << std::endl; // 
    writeSol(sol_path, bestSol.tour, bestSol.cost);
    writeTimeObj(save_path, argv[1], bestSol.cost, currTime);
    return 0;
}