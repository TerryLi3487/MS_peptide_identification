# MS_peptide_identification
Scoring function: 

To develop a scoring function, we need to calculate the similarity of mass spectrum data and the peptide data. 
I chose to discretize the mass array in the spectrum data and the fragments mass array for the peptide using a bin size of 1.0.
For each bin of mass, we take the max intensity of the spectrum data.
If a peptide fragment falls into that bin, then it has an intensity of 1 otherwise 0. 
Then, we will have two vectors of the same length, one for the intensity of spectrum data and one for peptide fragments. 
Then we divide each vector with its 2-norm and take the inner product of the two as the likelihood score. 


The number of identifications at 1% FDR:
Using the provided ups.fasta, with test1.mgf, I was able to get n=164 and with test2.mgf, I was able to get n=495.

Effort to optimize the search:
Instead of going through every spectrum and peptide pair, we compute the mass of the peptide the spectrum represents and only calculate the similarity if the difference between the mass of the peptide and the mass of the spectrum is less than 0.1Da. This way, we can save a lot of run time by not examining each pair. 
