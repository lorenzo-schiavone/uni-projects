#include "utilities.h"
#include <fstream>
#include <string>
#include <iostream>

// write solution
void writeSol(const char* filepathname, const std::vector<int>& tour, const double objval){
	std::ofstream File(filepathname);
    
    if (!File)
    {
      std::cout << "Something went wrong with opening the file!";
    }
    else
    { 
      File << objval << "\n";
      for ( uint i = 0; i < tour.size(); i++ ) {
      File << tour[i] << " ";
    }
        File.close();
    }
}

// write time and objval
void writeTimeObj(const char* filepathname, const char* instance, double objval, double time){
	std::ofstream File(filepathname, std::ios::app);
    
    if (!File)
    {
      std::cout << "Something went wrong with opening the file!";
    }
    else
    { 
      File << instance << "," << time << "," << objval << "\n" ;
      File.close();
    }
}

double getTimeLimit(const char* filepathname, const int n){

  std::ifstream file(filepathname);
    if (!file) {
        std::cerr << "Error: Could not open file.\n";
        return 1;
    }

    int size;
    double time_limit;

    while (file >> size >> time_limit) {
      if (size==n)
        return time_limit;
    }
    std::cerr << "Size not found.\n" << std::endl;
    return 0.0;
}

std::tuple<int, int, double, double, double, int> getGeneticParam(const char* filepathname, const int n){

  int size;
  int popSize = 0;
  int newPopsize = 0;
  double lsRate = 0.;
  double MAXmutationRate = 0.;
  double replace_factor = 0.;
  int MAXIT = 0;
  std::ifstream file(filepathname);
  if (!file) {
      std::cerr << "Error: Could not open file.\n";
      return std::make_tuple(popSize,newPopsize, lsRate,MAXmutationRate, replace_factor, MAXIT );
  }

  while (file >> size >> popSize >> newPopsize >> lsRate >> MAXmutationRate >> replace_factor >> MAXIT ) {
      if (size==n)
      return std::make_tuple(popSize, newPopsize, lsRate,MAXmutationRate, replace_factor, MAXIT );
  }
  std::cerr << "size not found" << std::endl;
  return std::make_tuple(popSize, newPopsize, lsRate, MAXmutationRate, replace_factor, MAXIT );
}