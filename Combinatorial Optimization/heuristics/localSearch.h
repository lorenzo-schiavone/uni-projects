#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include "TSPSolution.h"
#include "TSP.h"

struct TSPMove {
    int from;
    int to;
};

struct TSP3Move {
    int from;
    int mid;
    int to;
    int moveType;
};

bool localSearch(const TSPInstance& tsp, TSPSolution& sol, bool verbose=false);
bool localSearch2opt(const TSPInstance &tsp, TSPSolution &sol, bool verbose=false);


#endif // LOCAL_SEARCH_HPP