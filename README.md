<h3> Motif Finding </h3>
<p> 
Finds motifs within a set of DNA sequences. Created a benchmark program to generate synthetic data sets.
The benchmark program takes in 4 inputs: motif length, number of variable positions, sequence length, and sequence count. 
It then uses a random number generator to build sequences according to these parameters. 

Implemented two algorithms for motif finding: greedy and gibbs sampling. 

Modified the greedy algorithm to find
the closets motifs out of the first three sequences using the hamming distance. Build an initial
profile matrix out of these three motifs. Then goes through the rest of the DNA sequences and selects the best motif
in that sequence that maximizes profile matrix. 

Gibbs sampling:
Sampled sites from a probability distribution derived from PWM (position weight matrix). PWM is updated with
the frequencies of each character in the motifs. 

Evaluation Program: 
Compared number of overlapping sites with actual sites and running time.
</p>
