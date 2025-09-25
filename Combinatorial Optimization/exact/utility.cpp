#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdio.h>
#include <unordered_map>

using namespace std;

void read(const char* filename, vector<vector<double> >& C, int& N)
{
std::ifstream in(filename);
// read size
in >> N;
std::cout << "number of nodes n = " << N << std::endl;
// read costs
C.resize(N);
for (int i = 0; i < N; i++) {
	C[i].reserve(N);
	for (int j = 0; j < N; j++) {
	double c;
	in >> c;
	C[i].push_back(c);
	}
}
in.close();
}

void writeSol(const char* filepathname, const std::vector<int> tour, const double objval){
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

std::vector<int> reconstructTour(const std::vector<std::pair<int, int>>& edges) {
    std::unordered_map<int, int> nextNode; // dictionary
    
    for (const auto& edge : edges) {
        nextNode[edge.first] = edge.second;
    }

    std::vector<int> tour;
    int current = 0; 
    for (int i = 0; i < edges.size(); ++i) { 
        tour.push_back(current);
        current = nextNode[current]; 
    }
    tour.push_back(0); 
    return tour;
}

void writeTimeObj(const char* filepathname, const char* instance, double objval, double time){
	std::ofstream File(filepathname, std::ios::app);
    
    if (!File)
    {
      std::cout << "Something went wrong with opening the file!";
    }
    else
    { 
      File << instance << "," << objval << "," << time << "\n" ;
      File.close();
    }
}