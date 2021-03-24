import pandas as pd
from Bio import SeqIO
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description='Analyze a set of sample repeats for similarity and coverage statistics')
parser.add_argument('--multi_fasta', '-f', type=str, help='Multi-fasta to evaluate sequences from', required=True)
parser.add_argument('--snp_dists', '-s', type=str, help='CSV file of SNP differenes from snp-dists in molten format',
                    required=True)
parser.add_argument('--output_dir', '-o', type=str, help='output directory for single modified fasta files',
                    required=True)
parser.add_argument('--fasta_mixed', '-m', type=str, help='Multi-fasta without exclusively ACGT to count mixed positions',
                    required=True)

args = parser.parse_args()

snp_dists = pd.read_csv(args.snp_dists, header=None)
snp_dists.columns = ['sam_1', 'comparing', 'distance']

# filter if the sample 1 name is in the second column, and filter by snp distances
new_frame = snp_dists[[x[0] in x[1] for x in zip(snp_dists['sam_1'], snp_dists['comparing'])]]
identical_frame = new_frame.loc[(new_frame['sam_1'] != new_frame['comparing'])][snp_dists.distance == 0]
non_identical_frame = new_frame.loc[(new_frame['sam_1'] != new_frame['comparing'])][snp_dists.distance != 0]

fasta_sequences = SeqIO.parse(open(args.multi_fasta), 'fasta')

sam_n_counts = []
non_identical = []
for record in fasta_sequences:
# if the record is in the identical frame, include the N counts for comparison
# otherwise, add the ide to the lsit of non-identical
    if record.id in identical_frame.sam_1.unique() or record.id in identical_frame.comparing.unique():
        keys = ['sample_name', 'N_counts']
        values = [record.id, record.seq.count("N")]
        to_add = dict(zip(keys, values))
        sam_n_counts.append(to_add)
    else:
        non_identical.append(record.id)

n_counts_frame = pd.DataFrame(sam_n_counts)

os.chdir(args.output_dir)
if n_counts_frame.shape[0] != 0:
    n_counts_frame.sort_values(by=['sample_name'], ascending=True).to_csv("identical_sequences.csv", index=False)
if len(non_identical) != 0:
    with open("Non-identical sequences.txt", "w") as handle:
        for lines in non_identical:
            handle.write(lines+"\n")

print("Identical Sequences:")
if n_counts_frame.shape[0] != 0:
    print(n_counts_frame.sort_values(by=['sample_name'], ascending=True).to_string(index=False))
else:
    print("None")

print("Non-Identical sequences:")
if len(non_identical) != 0:
    for elem in non_identical:
        print(elem)
else:
	print("None")

# group the samples by their common name and find the first sample with the fesert N counts
n_counts_frame['standard_name'] = n_counts_frame['sample_name'].map(lambda x: x.rstrip('-v2|-v3|-v4|-v5|-v6'))
grouped = n_counts_frame.groupby(by=['standard_name'])
min_frame = n_counts_frame.loc[n_counts_frame.groupby('standard_name').N_counts.idxmin()][["sample_name", "N_counts"]]

fasta_mixed = SeqIO.parse(open(args.fasta_mixed), 'fasta')
nuc = ['A', 'T', 'C', 'G', 'N']
mixed_post = []
pos = {}
for record in fasta_mixed:
    mixed_pos = []
# if there are any mixed positions in the SNPs only fasta, append them to the data frame
    if record.id not in non_identical and record.id != "MN908947":
        count = 0
        for index, item in enumerate(record.seq):
            if item not in nuc:
                count += 1
                mixed_pos.append("{}:{}".format(item, int(index + 1)))

        keys = ['sample_name', 'mixed_counts']
        values = [str(record.id), count]
        to_add = dict(zip(keys, values))
        mixed_post.append(to_add)
        pos[str(record.id)] = ';'.join(n.strip("[]") for n in mixed_pos)


mixed_counts = pd.DataFrame(mixed_post)
positions_frame = pd.DataFrame(pos.items(), columns=['sample_name', 'mixed_positions'])
final_frame = mixed_counts.merge(positions_frame, on='sample_name', how='left').\
    merge(n_counts_frame, on='sample_name', how='left').drop(['standard_name'], axis=1)
# assign a designation for sample exclusion based on the N counts
final_frame['exclude_analysis'] = np.where(final_frame['sample_name'].isin(min_frame['sample_name']), "N", "Y")
if mixed_counts.shape[0] != 0:
    final_frame.sort_values(by=['sample_name']).to_csv("identical_sequences.csv", index=False)




