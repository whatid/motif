#include "benchmark.h" 
#include "motif.h"
#include <ctime>
#include <sstream>
#include <math.h>


int main(int argc, char*argv[])
{

	clock_t timer; 
		int ML, NM, SL, SC, num; 
		if (argc == 6) {
			ML = atoi(argv[1]);
			NM = atoi(argv[2]);
			SL = atoi(argv[3]);
			SC = atoi(argv[4]);		
			num = atoi(argv[5]); 
		} else {
			cout << "Error: invalid number of arguments" << endl;
			cout << "Arguments: motif length, number of variable positions, sequence length, sequence count" << endl; 
			return 0; 
		}

	srand(time(NULL)); 

	cout << "seed: " << time(NULL) << endl; 

	benchmark bench; 

	// default ML = 8, NM = 1, SL = 500, SC = 10
	bench.build(ML, NM, SL, SC, "data_set" + to_string(num));


		ifstream sequences; 
		sequences.open(("data_set" + to_string(num) + "/sequences.fa").c_str()); 
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
		motiflength.open(("data_set" + to_string(num) + "/motiflength.txt").c_str());
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

		ofstream runningtime; 
		runningtime.open("runningtime.txt", std::ios_base::app); 

		timer = clock(); 
		greedy.search(temp, parsed_ML); 
		runningtime << ((clock() - timer) / (double) CLOCKS_PER_SEC) * 1000 << endl; 


		ofstream output(("data_set" + to_string(num) + "/predictedsites.txt").c_str()); 

		vector<int> sites = greedy.predicted_sites; 
		for (int i = 0; i < sites.size(); i++)
		{
			output << sites[i] << endl; 
		}
		output.close(); 

		vector<string> seq = greedy.matrix; 

		ofstream PWM(("data_set" + to_string(num) + "/predictedmotif.txt").c_str());

		double PWM_matrix[4][parsed_ML];

		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < parsed_ML; y++)
			{
				PWM_matrix[x][y] = 0; 
			}
		}

		for (int i = 0; i < seq.size(); i++)
		{
			for (int j = 0; j < parsed_ML; j++)
			{
				if ((seq[i])[j] == 'A')
					PWM_matrix[0][j]++;
				else if ((seq[i])[j] == 'C')
					PWM_matrix[1][j]++; 
				else if ((seq[i])[j] == 'G')
					PWM_matrix[2][j]++; 
				else 
					PWM_matrix[3][j]++;  
			}
		}

		PWM << ">PMOTIF " << parsed_ML << endl; 
		for (int y = 0; y < parsed_ML; y++)
		{
			for (int x = 0; x < 4; x++)
			{
				PWM << PWM_matrix[x][y] << " "; 
				PWM_matrix[x][y] = (PWM_matrix[x][y] + (double(SC) / 1000)) / (double (SC) + (double(SC) / 250)); 
			}
			PWM << endl; 
		} 
		PWM << "<";

		PWM.close(); 

		double ** PWM_actual = bench.normalized_PWM; 

		double relative_entropy = 0; 

		for (int y = 0; y < parsed_ML; y++)
		{
			for (int x = 0; x < 4; x++)
			{
				relative_entropy += PWM_matrix[x][y] * log2(PWM_matrix[x][y] / PWM_actual[x][y]);  
			}
		} 

		ofstream entropy; 
		entropy.open("entropy.txt", std::ios_base::app); 

		entropy << relative_entropy << endl; 

		for (int i = 0; i < 4; i++)
		{
			delete [] PWM_actual[i]; 
		}
		delete[] PWM_actual; 

}
