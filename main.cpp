#include "benchmark.h" 
#include "motif.h"
#include <sstream>


int main(int argc, char*argv[])
{
	int ML, NM, SL, SC; 
		if (argc == 5) {
			ML = atoi(argv[1]);
			NM = atoi(argv[2]);
			SL = atoi(argv[3]);
			SC = atoi(argv[4]);		
		} else {
			cout << "Error: invalid number of arguments" << endl;
			cout << "Arguments: motif length, number of variable positions, sequence length, sequence count" << endl; 
			return 0; 
		}

	benchmark bench; 

	// default ML = 8, NM = 1, SL = 500, SC = 10
	bench.build(ML, NM, SL, SC, "data_set");	


	ifstream sequences; 
	sequences.open("data_set/sequences.fa"); 
	map<int, string> temp; 
	int count = -1; 

	for (string line; getline(sequences, line); )
	{
		if (line.find(">") == 0) {
			count++; 
			continue; 
		} else {
			temp[count].append(line);
		}
	}
	sequences.close(); 

	ifstream motiflength; 
	motiflength.open("data_set/motiflength.txt");
	int parsed_ML; 
	motiflength >> parsed_ML; 
	motiflength.close(); 

//  sanity check
/*
	for (int i = 0; i < 10; i++)
	{
		cout << "string" << i << endl; 
		cout << temp[i] << endl; 
		cout << "\n"; 
	}
*/
	motif greedy; 
	greedy.search(temp, parsed_ML); 
}
