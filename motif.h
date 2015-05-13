#ifndef MOTIF_H
#define MOTIF_H

#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>
#include <random>

using namespace std; 

class motif
{
	public: 
		vector<string> matrix; 
		vector<int>predicted_sites; 
		void search(map<int, string> DNA, int ML); 
		int score(string sequence, vector<string> prof_matrix);
};

#endif