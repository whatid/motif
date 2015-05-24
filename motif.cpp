#include "motif.h"

/* Gibbs Sampling */ 

void motif::search(map<int, string> DNA, int ML)
{
	random_device rd; 
	mt19937 gen(rd()); 

	// initialize predicted sites vector with random number generator. 
	for (int i = 0; i < DNA.size(); i++)
	{
		predicted_sites.push_back(rand() % (DNA[0].size() - ML));
	}

	// iterate until convergence. In reality, we let this program run for a very long time to achieve better results. 
	bool converge = false; 
	while (!converge)
	{
		string current (predicted_sites.begin(), predicted_sites.end()); 
		for (int index = 0; index < predicted_sites.size(); index++)
		{

		int rand_idx = index; //rand() % DNA.size(); 
		string rand_string (DNA[rand_idx]);  

		// zero out profile matrix 
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

		// create a profile matrix 
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
		
		// set predicted site based on probabilities 
		discrete_distribution<> dist(probabilities.begin(), probabilities.end());  
		predicted_sites[rand_idx] = dist(gen); 


		}

		// count the number of sites that have changed
		int num = 0; 
		for (int index = 0; index < predicted_sites.size(); index++)
		{
			if (current.at(index) != predicted_sites[index])
			{
				num++;  
			}
		}

		// this won't ever happen (no convergence). 
		if (num == 0)
			converge = true; 
		else 
			converge = false; 
	}
}

/*

#include "motif.h"


// Greedy Motif Search general outline 
// 1. Find 2 K-mers in sequences 1 and 2, form 2 x K alignment matrix and compute Score(S); pick the S with the highest score. 
// 2  Iteratively add 1 K-mer from each of the other (t-2) sequences.
// 3. At each of the following t-2 iterations, find a best K-mer in sequence i. 
 

void motif::search(map<int, string> DNA, int ML)
{

	for (int i = 0; i < DNA.size(); i++)
	{
		predicted_sites.push_back(0);
	}

	string string1 (DNA[0]);  
	string string2 (DNA[1]);
	string string3 (DNA[2]);

	string bestmotif1; 
	string bestmotif2;
	string bestmotif3;

	int lowest_distance = ML*3;

	for (int s1 = 0; s1 < string1.length()-ML; s1++)
	{
		string motif1 (&string1[s1], &string1[s1+ML]); 
	
		for (int s2 = 0; s2 < string2.length()-ML; s2++)
		{
			string motif2 (&string2[s2], &string2[s2+ML]); 
			for (int s3=0;s3<string3.length()-ML;s3++)
			{
				string motif3 (&string3[s3], &string3[s3+ML]);
				// calculate hamming distance
    			int distance = 0; 
				for (int i = 0; i < motif2.length(); i++)
				{
					if (motif1.at(i) != motif2.at(i))
						distance ++;
					if (motif2.at(i) != motif3.at(i))
						distance ++;
					if (motif1.at(i) != motif3.at(i))
						distance ++;
				}
				if (distance <= lowest_distance)
				{
					lowest_distance = distance; 
					bestmotif1 = motif1;
					bestmotif2 = motif2;
					bestmotif3 = motif3;

					predicted_sites[0] = s1; 
					predicted_sites[1] = s2; 
					predicted_sites[2] = s3;
				}
			}      	
		}
	}
		matrix.push_back(bestmotif1); 
		matrix.push_back(bestmotif2); 
		matrix.push_back(bestmotif3);


	for (int i = 3; i < DNA.size(); i++)
	{
	//	cout << i << endl; 
		int best_score = 0;
		string best_motif_i;
		string string_i (DNA[i]);
		for (int s_i = 0; s_i < string_i.length()-ML; s_i ++)
		{
			
			string motif_i (&string_i[s_i], &string_i[s_i + ML]);
			matrix.push_back(motif_i); 
			// calc score
			int temp_score = score(motif_i, matrix);

			if (temp_score > best_score)
			{
				best_score = temp_score;
				best_motif_i = motif_i;  
				predicted_sites[i] = s_i; 
			}

			matrix.pop_back(); 
		}
		matrix.push_back(best_motif_i);
	}

}

int motif::score(string sequence, vector<string> prof_matrix)
{

	int width = prof_matrix[0].length(); 
	int length = 4; 

	int score_matrix[width][length];

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < length; y++)
		{
			score_matrix[x][y] = 0; 
		}
	}

	for (int i = 0; i < prof_matrix.size(); i++)
	{
		int A = 0, C = 0, G = 0, T = 0; 
		for (int j = 0; j < prof_matrix[0].length(); j++)
		{
			if (prof_matrix[i].at(j) == 'A')
				score_matrix[j][1] += 1; 
			else if (prof_matrix[i].at(j) == 'C')
				score_matrix[j][2] += 1; 
			else if (prof_matrix[i].at(j) == 'G')
				score_matrix[j][3] += 1; 
			else 
				score_matrix[j][4] += 1; 
		}
	}

	int final_score = 0; 

	for (int x = 0; x < width; x++)
	{
		final_score += max(max(score_matrix[x][1], score_matrix[x][2]) , max(score_matrix[x][3], score_matrix[x][4])); 
	}

	return final_score; 
}
*/


















