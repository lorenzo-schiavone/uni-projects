#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include "cpxmacro.h"
#include <random>
#include <string>
#include <sys/time.h>
#include <stdio.h>
#include "utility.cpp"

using namespace std; 
 
// error status and messagge buffer (from cpxmacro.h)
int status;
char errmsg[BUF_SIZE];
const int NAME_SIZE = 512;
char name[NAME_SIZE]; // for varaible names
double zero = 0.0;
double one = 1.0;

// data
int N; //number of nodes
vector<vector<double> > C; // cost matrix
double infinite; // infinite value (an upper bound on the value of any feasible solution)
vector<vector<int> > map_x;	// x_ij ---> map_x[i][j]	
vector<vector<int> > map_y;	// y_ij ---> map_y[i][j]	
std::vector<std::pair<int, int> > rev_map; // indice -> i,j

int idx_start_y; // index to know when y variables start
double epsilon = 1e-8;

void setupLP(CEnv env, Prob lp)
{	
    int current_var_position = 0; 
	// map_x N by N initialized by -1
	map_x.resize(N, vector<int>(N,-1)); // with cpp11 works, for older version use the below
	// map_x.resize(N);
	// for ( int i = 0 ; i < N ; ++i ) {
	// 	map_x[i].resize(N, -1);
	// 	// for ( int j = 0 ; j < N ; ++j ) {
	// 	// 	map_x[i][j] = -1;
	// 	// }
	// } 
	for (int i = 0; i < N; i++)
	{
		for (int j = 1; j < N; j++) // j starts from 1 to avoid creating more variable than needed: C_ii = 0 always
		{
			if ( C[i][j]<epsilon){continue;} //EXT1 -- no od_cost_max yet -> C[i][j] > od_cost_max 
			char xtype = 'I'; // integer
			double lb = 0.0;
			double ub = CPX_INFBOUND; 
			snprintf(name, NAME_SIZE, "x_%d_%d", i, j);
			char* xname = (char*)(&name[0]);
			
			// crete variable one by one
			CHECKED_CPX_CALL( CPXnewcols, env, lp, 1   , &(zero), &lb, &ub, &xtype, &xname ); // zero as it doesn't appear in the obj function
			/// status =      CPXnewcols (env, lp, ccnt, obj      , lb  , ub, xctype, colname);
			
			map_x[i][j] = current_var_position ++;
			rev_map.push_back(std::make_pair(i,j));

		}
	}
	// std::cout << "x created" << std::endl;

	idx_start_y = rev_map.size();
	// map_y N by N initialized by -1
	map_y.resize(N, vector<int>(N,-1)); // with cpp11 works, for older version use the below
	// map_y.resize(N);
	// for ( int i = 0 ; i < N ; ++i ) {
	// 	map_y[i].resize(N);
	// 	for ( int j = 0 ; j < N ; ++j ) {
	// 		map_y[i][j] = -1;
	// 	}
	// } 

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if ( C[i][j]<epsilon) {continue;} //C[i][j] > od_cost_max
			char ytype = 'B'; // binary type
			double lb = 0.0;
			double ub = 1.0; 
			snprintf(name, NAME_SIZE, "y_%d_%d", i, j);
			char* yname = (char*)(&name[0]);
			// crete variable one by one
			CHECKED_CPX_CALL( CPXnewcols, env, lp, 1  , &(C[i][j]), &lb, &ub, &ytype, &yname );
			/// status =      CPXnewcols (env, lp, ccnt, obj      , lb  , ub, xctype, colname);

			map_y[i][j] = current_var_position ++;
			rev_map.push_back(std::make_pair(i,j));
		}
	}
	// std::cout << "y created" << std::endl;

	// CONSTRAINTS: 

	// constr1
	for (int k = 1; k < N; k++)
	{
		std::vector<int> idx;
		std::vector<double> coef;
		char sense = 'E';
		for (int i = 0; i < N; i++)
		{
			if ( map_x[i][k] < 0 ) {continue;}  				///EXT1
			
			idx.push_back(map_x[i][k]); 
			coef.push_back(1.0); 
		}
		for (int j = 0; j < N; j++)
		{
			if ( map_x[k][j] < 0 ) {continue;}  				///EXT1
			
			idx.push_back(map_x[k][j]); 
			coef.push_back(-1.0); 
		}
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , idx.size(), &(one), &sense, &matbeg, &idx[0], &coef[0], NULL      , NULL      );
	}
	// std::cout << "constraint 1 created" << std::endl;

	// constr2
	for (int i = 0; i < N; i++){
		std::vector<int> idx;
		std::vector<double> coef;
		char sense = 'E';
		for (int j=0; j<N; j++)
		{
			if ( map_y[i][j] < 0 ) continue;  				///EXT1
			
			idx.push_back(map_y[i][j]);
			coef.push_back(1.0); 
		}
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , idx.size(), &(one), &(sense), &matbeg, &idx[0], &coef[0], NULL      , NULL      );
	} 
	// std::cout << "constraint 2 created" << std::endl;

	// constr3
	for (int j = 0; j < N; j++){
		std::vector<int> idx;
		std::vector<double> coef;
		char sense = 'E';
		for (int i=0; i<N; i++)
		{
			if ( map_y[i][j] < 0 ) continue;  				///EXT1
			
			idx.push_back(map_y[i][j]);
			coef.push_back(1.0); 
		}
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , idx.size(), &(one), &(sense), &matbeg, &idx[0], &coef[0], NULL      , NULL      );
	}
	// std::cout << "constraint 3 created" << std::endl;

	// constr4
	for (int i = 0; i < N; i++){
		char sense = 'L';
		for (int j=0; j<N; j++)
		{
			std::vector<int> idx;
			std::vector<double> coef;
			if ( (map_y[i][j] < 0) || (map_x[i][j] < 0)) continue;  				///EXT1
			
			idx.push_back(map_y[i][j]);
			coef.push_back(1.0-static_cast<double>(N));

			idx.push_back(map_x[i][j]);
			coef.push_back(1.0); 
			int matbeg = 0;
			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , idx.size(), &(zero), &(sense), &matbeg, &idx[0], &coef[0], NULL      , NULL      );
		}
		
	} 

	// to check if everything is right
	// CHECKED_CPX_CALL( CPXwriteprob, env, lp, "tsp.lp", NULL );
}

int main (int argc, char const *argv[])
{	
	try
		{
		// Base paths
		std::string base_path1 = "../data/cost_matrix/cost";
		std::string base_path2 = "../solutions/exact/sol";

		// File extensions
		std::string extension1 = ".dat";
		std::string extension2 = ".txt";

		// Construct the full paths
		std::string path1 = base_path1 + argv[1] + extension1;
		std::string path2 = base_path2 + argv[1] + extension2; 

		const char* costs_path = path1.c_str();
    	const char* sol_path = path2.c_str();

		std::string base_path3 = "../results/exact_summary.csv";
		const char* summ_path = base_path3.c_str();

		read(costs_path, C, N);

		DECL_ENV( env ); 
		DECL_PROB( env, lp ); 
		
		std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
		///// setup MIP
		setupLP(env, lp);
		// std::cout << "setup completed" << std::endl;

	// 	///// optimize
		CHECKED_CPX_CALL( CPXmipopt, env, lp ); 
		std::chrono::high_resolution_clock::time_point finish= std::chrono::high_resolution_clock::now();
		std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

		// print objective function value
		std::cout << "Optimization Completed" << std::endl;
		std::cout << "in " << duration.count()/1000.0 << " seconds (CPU time)\n";
		
		double objval;
		CHECKED_CPX_CALL( CPXgetobjval, env, lp, &objval ); // get the objective value
		std::cout << "Objval: " << objval << std::endl;
		
		// ALTERNATIVE TO TAKE ONLY THE YS DIRECTLY
		std::vector<double> varVals; // vector of doubles where we store what we want to show
		int n = CPXgetnumcols(env,lp); // number of variabels
		varVals.resize(n-idx_start_y); // take all
		std::vector<std::pair<int, int> > result;
		CHECKED_CPX_CALL( CPXgetx ,env, lp, &varVals[0], idx_start_y, n-1);
		for (int p=0; p<n-idx_start_y; p++){ // only the ys. if 0 no connection, if > 0 there is the edge in the tour
			if (varVals[p] <= epsilon ) {continue;}
			result.push_back(rev_map[p+idx_start_y]);
		}
		
		std::vector<int> tour = reconstructTour(result);
		writeSol(sol_path, tour, objval);
		writeTimeObj(summ_path, (const char* ) argv[1], objval, (double) duration.count()/1000);
			
		///// free memory
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}

