#include "motif.h"

/* Greedy Motif Search general outline 
1. Find 2 K-mers in sequences 1 and 2, form 2 x K alignment matrix and compute Score(S); pick the S with the highest score. 
2  Iteratively add 1 K-mer from each of the other (t-2) sequences.
3. At each of the following t-2 iterations, find a best K-mer in sequence i. 
*/ 

void motif::search(map<int, string> DNA, int ML)
{
	random_device rd; 
	mt19937 gen(rd()); 
	for (int i = 0; i < DNA.size(); i++)
	{
		predicted_sites.push_back(rand() % (DNA[0].size() - ML));
	}
	bool converge = false; 
	while (!converge)
	{
		string noob (predicted_sites.begin(), predicted_sites.end()); 
		for (int index = 0; index < predicted_sites.size(); index++)
		{

		int rand_idx = index; //rand() % DNA.size(); 
		string rand_string (DNA[rand_idx]);  

		double profile[4][ML];

		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < ML; y++)
			{
				profile[x][y] = 0.0	; 
			}
		}

		for (int i = 0; i < DNA.size(); i++)
		{
			if (i == rand_idx)
				continue; 
			else 
			{
				string lmer (&(DNA[i])[predicted_sites[i]], &(DNA[i])[predicted_sites[i] + ML]); 

				//cout << lmer << endl; 
				for (int j = 0; j < ML; j++)
				{
					if (lmer.at(j) == 'A')
						profile[0][j]++;
					else if (lmer.at(j) == 'C')
						profile[1][j]++; 
					else if (lmer.at(j) == 'G')
						profile[2][j]++; 
					else 
						profile[3][j]++;  
				}
			}
		}

		for (int y = 0; y < ML; y++)
		{
			for (int x = 0; x < 4; x++)
			{
				profile[x][y] = double(profile[x][y]) / double(9); 
			}
		} 

		map <int, double> distribution; 
		double lowest = INFINITY; 
		for (int i = 0; i < rand_string.size()-ML; i++)
		{	
			double result = 1; 
			string motif_i(&rand_string[i], &rand_string[i+ML]); 
			for (int j = 0; j < motif_i.length(); j++)
			{
				if (motif_i.at(j) == 'A')
					result *= profile[0][j]; 
				else if (motif_i.at(j) == 'C')
					result *= profile[1][j]; 
				else if (motif_i.at(j) == 'G')
					result *= profile[2][j]; 
				else 
					result *= profile[3][j]; 
			}
			distribution[i] = result; 

			if (result != 0 && result <= lowest)
				lowest = result; 
		}
		//cout << lowest << endl; 
		double sum = 0.0; 
		for (int i = 0; i < distribution.size(); i++)
		{
			distribution[i] = distribution[i] / lowest; 
			sum += distribution[i]; 
			//cout << sum << endl; 
		}

		vector<double> probabilities; 
		for (int i = 0; i < distribution.size(); i++)
		{
			distribution[i] = distribution[i] / sum; 
			probabilities.push_back(distribution[i]); 
		}
		
		discrete_distribution<> dist(probabilities.begin(), probabilities.end());  
		predicted_sites[rand_idx] = dist(gen); 


		}

		bool breh = false; 
		for (int index = 0; index < predicted_sites.size(); index++)
		{
			if (noob.at(index) != predicted_sites[index])
			{
				breh = true; 
				break; 
			}
		}

		if (breh)
			converge = false; 
		else 
			converge = true; 
	}
}









