#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <vector>
#include <fstream>
#include <string>

using namespace std; 
class benchmark
{
	public: 
		void build(int ML, int NM, int SL, int SC, string folder); 
		map <int, vector<char> > randSequence; 
		vector<char> motif; 
		map <int, vector<char> > bindingSites; 
		vector<int> siteLocation; 
};
#endif 