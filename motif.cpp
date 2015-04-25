#include "motif.h"

/* Greedy Motif Search general outline 
1. Find 2 K-mers in sequences 1 and 2, form 2 x K alignment matrix and compute Score(S); pick the S with the highest score. 
2  Iteratively add 1 K-mer from each of the other (t-2) sequences.
3. At each of the following t-2 iterations, find a best K-mer in sequence i. 
*/ 

void motif::search(map<int, string> DNA, int ML)
{
	string string1 (DNA[0]);  
	string string2 (DNA[1]); 

	string bestmotif1; 
	string bestmotif2; 

	int lowest_distance = ML; 

	for (int s1 = 0; s1 < string1.length()-ML; s1++)
	{
		string motif1 (&string1[s1], &string1[s1+ML]); 
	
		for (int s2 = 0; s2 < string2.length()-ML; s2++)
		{
			string motif2 (&string2[s2], &string2[s2+ML]); 

			// calculate hamming distance
			int distance = 0; 
			for (int i = 0; i < motif2.length(); i++)
			{
				if (motif1.at(i) != motif2.at(i))
					distance ++; 
			}
			if (distance < lowest_distance)
			{
				lowest_distance = distance; 
				bestmotif1 = motif1; 
				bestmotif2 = motif2; 
			}
		}
	}

	// form K x 4 seed matrix
	vector<string> matrix; 
	matrix.push_back(bestmotif1); 
	matrix.push_back(bestmotif2); 

	for (int i = 2; i < DNA.size(); i++)
	{
		int best_score = 0; 
		string best_motif_i = NULL; 
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

	fill(&score_matrix[0][0], &score_matrix[0][0] + sizeof(score_matrix), 0); 

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









