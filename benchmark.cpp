#include "benchmark.h"



void benchmark::build(int ML, int NM, int SL, int SC, string folder)
{
	
	// step 1

	srand(time(NULL));  

	char nucleotide [4] = {'A', 'C', 'T', 'G'};

	//step 2

	for (int i = 0; i < SC; i++)
	{
		for (int j = 0; j < SL; j++)
		{
			randSequence[i].push_back(nucleotide[rand() % 4]);
		}
	}


	// step 3
	

	for (int i = 0; i < ML; i++)
	{
		motif.push_back(nucleotide[rand() % 4]);
	}

	for (int i = 0; i < NM; i++)
	{
		int idx = rand() % motif.size();

		// mark as variable
		motif.at(idx) = 'V';  
	}

	//step 4

	

	for (int i = 0; i < SC; i++)
	{
		for (int j = 0; j < ML; j++)
		{
			if (motif[j] == 'V')
				bindingSites[i].push_back(nucleotide[rand() % 4]); 
			else 
				bindingSites[i].push_back(motif[j]); 
		}
	}

	//step 5
	// put motif in sequence entirely 



	for (int i = 0; i < SC; i++)
	{
		int idx = rand() % (SL - ML + 1);
		siteLocation.push_back(idx);

		for (int j = idx; j < ML + idx; j++)
		{
			randSequence[i].at(j) = bindingSites[i].at(j - idx);
		}
		
	}

	//create folder
	system(("mkdir " + folder).c_str());

	//step 6
	ofstream seq; 
	seq.open((folder + "/sequences.fa").c_str());

	for (int i = 0; i < SC; i++)
	{
		seq << "> sequence " << i << endl;
		copy(randSequence[i].begin(), randSequence[i].end(), ostream_iterator<char>(seq, ""));  
		seq << endl; 
	}

	seq.close(); 


	// step 7
	ofstream sites; 
	sites.open((folder + "/sites.txt").c_str());
	for (int i = 0; i < siteLocation.size(); i++)
	{
		sites << siteLocation[i] << endl;
	}

	sites.close(); 


	// step 8 

	ofstream motifFile; 
	motifFile.open((folder + "/motif.txt").c_str());

	motifFile << "MOTIF " << ML << " "; 

	for (int i = 0; i < ML; i++)
	{
		if (motif[i] == 'V')
			motifFile << "*"; 
		else 
			motifFile << motif[i]; 
	}

	motifFile.close(); 


	//step 9 
	ofstream motifLength; 
	motifLength.open((folder + "/motiflength.txt").c_str());
	motifLength << ML; 
	motifLength.close(); 


}