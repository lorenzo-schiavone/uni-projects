#ifndef TSP_INSTANCE_HPP
#define TSP_INSTANCE_HPP

#include <vector>
#include <string>

struct TSPInstance {
    std::vector<std::vector<double>> costMatrix;
    int numCities;
    double infinite;

    void read(const char* filename);
};
#endif // TSP_INSTANCE_HPP