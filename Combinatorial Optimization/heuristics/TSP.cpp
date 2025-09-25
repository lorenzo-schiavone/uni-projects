#include "TSP.h"
#include <fstream>
#include <iostream>

void TSPInstance::read(const char* filename) {
    std::ifstream in(filename);
    in >> numCities;
    costMatrix.resize(numCities);
    for (int i = 0; i < numCities; i++) {
        costMatrix[i].reserve(numCities);
        for (int j = 0; j < numCities; j++) {
            double c;
            in >> c;
            costMatrix[i].push_back(c);
        }
    }
    in.close();

    infinite = 0;
    for (int i = 0; i < numCities; i++) {
        for (int j = 0; j < numCities; j++) {
            infinite += costMatrix[i][j];
        }
    }
    infinite *= 2;
}