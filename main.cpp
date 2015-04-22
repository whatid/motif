#include "benchmark.h" 
#include "motif.h"


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

	bench.build(ML, NM, SL, SC, "data_set");	
}
