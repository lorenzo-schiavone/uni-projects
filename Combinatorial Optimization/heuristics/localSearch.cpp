#include "localSearch.h"
#include "utilities.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <numeric>

void apply2optMove(TSPSolution& tspSol, const TSPMove& move) {
    TSPSolution tmpSol = tspSol;
    for (int i = move.from; i <= move.to; ++i) {
        tspSol.tour[i] = tmpSol.tour[move.to - (i - move.from)];
    }
}

void apply3optMove(TSPSolution& tspSol, const TSP3Move& move)
  {
    TSPSolution tmpSol(tspSol);

    int a = move.from; int b = move.mid; int c = move.to; int moveType = move.moveType;

    int index = 0;
    // first segment (0 to a)
    for (int i = 0; i <= a; ++i) {
        tspSol.tour[index++] = tmpSol.tour[i];
    }
    
    switch (moveType)
    {
        case 0:

        // From b to a+1, reversed
        for (int i = b; i >= a + 1; --i) {
            tspSol.tour[index++] = tmpSol.tour[i];  
        }
        // From c to b+1, reversed)
        for (int i = c; i >= b + 1; --i) {
            tspSol.tour[index++] = tmpSol.tour[i];
        }
        break;

        case 1:

        // From b+1 to c, in original order
        for (int i = b + 1; i <= c; ++i) {
            tspSol.tour[index++] = tmpSol.tour[i];
        }
        // From a+1 to b, in original order
        for (int i = a + 1; i <= b; ++i) {
            tspSol.tour[index++] = tmpSol.tour[i];
        }
        break;

        case 2:
        // From b+1 to c, in original order
        for (int i = b + 1; i <= c; ++i) {
            tspSol.tour[index++] = tmpSol.tour[i];
        }
        // From b to a+1, reversed
        for (int i = b; i >= a + 1; --i) {
            tspSol.tour[index++] = tmpSol.tour[i];  
        }
        break;

        case 3:
        // From c to b+1, reversed)
        for (int i = c; i >= b + 1; --i) {
            tspSol.tour[index++] = tmpSol.tour[i];
        }

        // From a+1 to b, in original order
        for (int i = a + 1; i <= b; ++i) {
            tspSol.tour[index++] = tmpSol.tour[i];
        }
        break;
        
        default:
        throw std::invalid_argument("Invalid moveType in apply3optMove");
        break;
    }
    for (int i = c + 1; i < tmpSol.tour.size(); ++i) {
        tspSol.tour[index++] = tmpSol.tour[i];
    }
  } 

// 2-opt neighbor search
double findBestNeighbor(const TSPInstance& tsp, const TSPSolution& currSol, TSPMove& move) {
    double bestCostVariation = tsp.infinite;
    for (int a = 1; a < currSol.tour.size() - 2; a++) {
        int h = currSol.tour[a - 1];
        int i = currSol.tour[a];
        for (int b = a + 1; b < currSol.tour.size() - 1; b++) {
            int j = currSol.tour[b];
            int l = currSol.tour[b + 1];
            double neighCostVariation = - tsp.costMatrix[h][i] - tsp.costMatrix[j][l] + tsp.costMatrix[h][j] + tsp.costMatrix[i][l];
            if (neighCostVariation < bestCostVariation) {
                bestCostVariation = neighCostVariation;
                move.from = a;
                move.to   = b;
            }
        }
    }
    return bestCostVariation;
}

double findBestNeighbor3opt(const TSPInstance& tsp, const TSPSolution& currSol, 
    TSP3Move& bestMove, double early_break) {
    double bestCostVariation = tsp.infinite;
    int n = currSol.tour.size();
    std::vector<double> increments(4);
    
    for (int a = 1; a < n - 3; a++) {
        int h = currSol.tour[a];
        int i = currSol.tour[a + 1];
        for (int b = a + 2; b < currSol.tour.size() - 2; b++) {
            int j = currSol.tour[b];
            int l = currSol.tour[b + 1];
            for (int c = b + 2; c < currSol.tour.size() - 1; c++) {
                int k = currSol.tour[c];
                int m = currSol.tour[c + 1];
                double moveCost = -tsp.costMatrix[h][i] - tsp.costMatrix[j][l] - tsp.costMatrix[k][m];
                increments[0] = tsp.costMatrix[h][j] + tsp.costMatrix[i][k] + tsp.costMatrix[l][m];
                increments[1] = tsp.costMatrix[h][l] + tsp.costMatrix[k][i] + tsp.costMatrix[j][m];
                increments[2] = tsp.costMatrix[h][l] + tsp.costMatrix[k][j] + tsp.costMatrix[i][m];
                increments[3] = tsp.costMatrix[h][k] + tsp.costMatrix[i][l] + tsp.costMatrix[j][m];
                int moveType = arg_min(increments);
                double totalVariation = moveCost + increments[moveType];
                if (totalVariation < bestCostVariation) {
                    bestCostVariation = totalVariation;
                    bestMove.from = a;
                    bestMove.mid  = b;
                    bestMove.to   = c;
                    bestMove.moveType = moveType;
                }
                if (bestCostVariation < early_break){
                    break; // first improvement
                }
            }
        }
    }
    return bestCostVariation;
}

// Local Search
bool localSearch(const TSPInstance &tsp, TSPSolution &sol, bool verbose) {
    if(verbose){
        std::cout << "\n";
    }
    try {
        bool stop = false;
        int iter = 0;
        double currValue = sol.cost; 
        double bestValue = currValue;
        TSPMove move;
        TSP3Move move3;
        double lastImprovVal= -1e12;
        
        while (!stop) {
            iter++;
            double improvement2opt = findBestNeighbor(tsp, sol, move);
            if (currValue + improvement2opt < bestValue - 1e-14) {
                currValue = currValue + improvement2opt;
                lastImprovVal = improvement2opt;
                // std::cout << lastImprovVal << std::endl;
                if (verbose){
                    std::cout <<"2-opt move:\n";
                    std::cout << " (" << ++iter << "ls) value " << currValue << " (from " << bestValue << ")\n";
                    std::cout << " move: " << move.from << " , " << move.to << std::endl;
                }
                bestValue = currValue;
                apply2optMove(sol, move);
                
            } else {
                double improvement3opt = findBestNeighbor3opt(tsp, sol, move3, .75*lastImprovVal); // .8*lastImprovVal -1e10
                // EARLY BREAK
                if (currValue + improvement3opt < bestValue - 1e-14) {
                    lastImprovVal = improvement3opt;

                    // std::cout << lastImprovVal << std::endl;    
                    currValue = currValue + improvement3opt;
                    if (verbose){
                        std::cout <<"3-opt move:\n";
                        std::cout << " (" << ++iter << "ls) value " << currValue << " (from " << bestValue << ")\n";
                        std::cout << " move" << move3.moveType << ": " << move3.from << " , " << move3.mid << " , " << move3.to << std::endl;
                    }
                    bestValue = currValue;
                    apply3optMove(sol, move3);
                }
                else {
                    stop = true;
                    // std::cout << "iter: " << iter << std::endl;
                }
            }
        }
        sol.cost = bestValue;
    } catch (std::exception &e) {
        std::cout << ">>>EXCEPTION in localSearch: " << e.what() << std::endl;
        return false;
    }

    return true;
}

bool localSearch2opt(const TSPInstance &tsp, TSPSolution &sol, bool verbose) {
    if(verbose){
        std::cout << "\n";
    }
    try {
        bool stop = false;
        int iter = 0;
        double currValue = sol.cost; 
        double bestValue = currValue;
        TSPMove move;
        while (!stop) {
            iter++;
            double improvement2opt = findBestNeighbor(tsp, sol, move);
            if (currValue + improvement2opt < bestValue - 1e-14) {
                currValue = currValue + improvement2opt;
                if (verbose){
                    std::cout <<"2-opt move:\n";
                    std::cout << " (" << ++iter << "ls) value " << currValue << " (from " << bestValue << ")\n";
                    std::cout << " move: " << move.from << " , " << move.to << std::endl;
                }
                bestValue = currValue;
                apply2optMove(sol, move);
                
            } else {
                stop = true;

            }
        }
        sol.cost = bestValue;
    } catch (std::exception &e) {
        std::cout << ">>>EXCEPTION in localSearch: " << e.what() << std::endl;
        return false;
    }
    return true;
}