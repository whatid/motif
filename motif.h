#ifndef MOTIF_H
#define MOTIF_H

#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <string>

using namespace std; 

class motif
{
	public: 
		void search(map<int, string> DNA, int ML); 
		int score(string sequence, vector<string> prof_matrix);
};

#endif