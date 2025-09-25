#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <vector>
#include <string>
#include <random>
#include "TSPSolution.h"

void writeSol(const char* filepathname, const std::vector<int>& tour, double objval);
void writeTimeObj(const char* filepathname, const char* instance, double objval, double time);
double getTimeLimit(const char* filepathname, const int n);
std::tuple<int, int, double, double, double, int> getGeneticParam(const char* filepathname, const int n);

// https://jclay.github.io/dev-journal/simple_cpp_argmax_argmin.html
template <typename T, typename A>
int arg_min(const std::vector<T, A>& vec) {
    return static_cast<int>(std::distance(vec.begin(), std::min_element(vec.begin(), vec.end())));
}

#endif // UTILITIES_HPP