# Script for ncov sample repeat analysis

Scripts for the comparison of repeat samples, which are internally designated in the following way: 
**sample_name-v#**, where # is the number of times the sample has repeated sequencing

## Identification of repeat samples

The script <em>identify_repeat_samples.py</em> will parse a directory containing output files from the ncov19 pipeline and identify any samples with the "-v#" suffix as repeats. It can then egnerate a .txt file of all repeat names and the
original samples. So for example, if sample_1-v3 is in the directory, the list will also populated with sample_1 and sample_1-v2. 


## Repeat Analysis 

Repeat analysis compares the common SNPs among samples with similar names to classify similarity. In the following list:

- sample_1
- sample_1-v2
- sample_2
- sample_2-v2

The script will consider the pair of sample_1 for comparison, and sample_2.

The following software dependencies are required for the <em>analyze_repeats.py</em> script:


- faSomeRecords
- MAFFT
- snp-sites>=2.5.1
- snp-dists>=0.7.0

### Similarity Logic

The script will categorize the pairs of samples in the following way: 

- If the SNPS at mutually covered sites are the same (disregarding N's), they will be classified as identical.
Otherwise, they will be classified as non-identical. 
- If two samples that are not part of a pair (i.e. sample_1-v2 and sample_2) have identical SNPs as above, they
will be categorized as identical, but different names. 

- For identical samples: 
	- Of at least two samples in the repeat group have genome completeness of 85% of greater, then the
	sample with fewer mixed positions will be kept for analysis. Otherwise, the sample with the greatest
	genome completeness will be kept for analysis (the number of mixed sites is used to break a tie)
	- of the number of mixed sites is greater than 2 for a sample, it is flagged with this detail

- For non-identical samples: 
	- Both samples are flagged for review for a possible sample mix-up if the number of mixed sites is greater
	than 2
	

The script also considers the overall coverage of a sample (measured in the concensus FASTA with an N count)
as well as the number of mixed SNP sites measured by snp-sites. If a SNP is not called an A, T, C, G, it will
be classified as a mixed site. 

## Outputs

The script has the following possible outputs: 

- all_repeats.csv: For all samples involved in repeat analysis, containing mixed sites, number of
mixed sites, and N count
- all_repeats_w_metadata.csv: For all samples involved in repeat analysis, containing mixed sites, number of
mixed sites, and N count. If the metadata file is provided from the database, the results are merged together
- non-identical_sequences.txt: If the samples have similar names but non-identical SNPs, their names are listred here.
This may be indicative of a sample swap or contamination. 
- identical_different_names.csv: For samples that are different names but identical SNPs, the output from the snp-dists
is exported here. 



