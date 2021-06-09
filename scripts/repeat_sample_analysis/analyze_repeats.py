import pandas as pd
from Bio import SeqIO
import argparse
import numpy as np
import warnings

warnings.filterwarnings("ignore")


def get_grouping(data_frame):
    data_frame['standard_name'] = data_frame['sample_name'].map(lambda x: x.split('-v', 1)[0])
    min_frame = data_frame.loc[data_frame.groupby('standard_name').N_counts.idxmin()][["sample_name", "N_counts"]]
    return data_frame, min_frame


def main():
    parser = argparse.ArgumentParser(description='Analyze a set of sample repeats for similarity and coverage'
                                                 'statistics')
    parser.add_argument('--multi_fasta', '-f', type=str, help='Multi-fasta to evaluate sequences from', required=True)
    parser.add_argument('--snp_dists', '-s', type=str, help='CSV file of SNP differences from snp-dists in molten'
                                                            'format', required=True)
    parser.add_argument('--output_dir', '-o', type=str, help='output directory for single modified fasta files',
                        required=True)
    parser.add_argument('--fasta_mixed', '-m', type=str, help='Multi-fasta without exclusively ACGT to count mixed'
                                                              'positions', required=True)
    parser.add_argument('--metadata', '-d', type=str, help='Optional metadata frame with additional sample information',
                        required=False)

    args = parser.parse_args()

    snp_dists = pd.read_csv(args.snp_dists, header=None)
    snp_dists.columns = ['sam_1', 'comparing', 'distance']

    # filter if the sample 1 name is in the second column, and filter by snp distances
    new_frame = snp_dists[[x[0] in x[1] for x in zip(snp_dists['sam_1'], snp_dists['comparing'])]]
    identical_frame = new_frame.loc[(new_frame['sam_1'] != new_frame['comparing'])][snp_dists.distance == 0]
    # non_identical_frame = new_frame.loc[(new_frame['sam_1'] != new_frame['comparing'])][snp_dists.distance != 0]

    diff_name_identical = snp_dists[~snp_dists.sam_1.isin(identical_frame.sam_1)]
    diff_name_identical = diff_name_identical[~diff_name_identical.sam_1.isin(identical_frame.comparing)]
    diff_name_identical = diff_name_identical[snp_dists.apply(lambda x: x.sam_1 not in
                                                                        x.comparing, axis=1)][snp_dists.distance == 0]
    diff_name_identical = diff_name_identical[~diff_name_identical['sam_1'].isin(['MN908947'])]

    if len(diff_name_identical) != 0:
        diff_name_identical.to_csv("different_names_identical.csv", index=False)

    fasta_sequences = SeqIO.parse(open(args.multi_fasta), 'fasta')

    sam_n_counts = []
    non_identical = []
    for record in fasta_sequences:
        if record.id in identical_frame.sam_1.unique() or record.id in identical_frame.comparing.unique():
            keys = ['sample_name', 'N_counts']
            values = [record.id, record.seq.count("N")]
            to_add = dict(zip(keys, values))
            sam_n_counts.append(to_add)
        else:
            non_identical.append(record.id)

    n_counts_frame = pd.DataFrame(sam_n_counts)
    non_identical.sort()

    if len(non_identical) != 0:
        with open("Non-identical sequences.txt", "w") as handle:
            for lines in non_identical:
                handle.write(lines+"\n")

    print("Identical Sequences:")
    if n_counts_frame.shape[0] != 0:
        print(n_counts_frame.sort_values(by=['sample_name'], ascending=True))
    else:
        print("None")

    print("Non-Identical sequences:")
    if len(non_identical) != 0:
        for elem in non_identical:
            print(elem)
    else:
        print("None")

    fasta_mixed = SeqIO.parse(open(args.fasta_mixed), 'fasta')
    # prevent deletions from being included with -
    nuc = ['A', 'T', 'C', 'G', 'N', '-']
    mixed_post = []
    pos = {}
    for record in fasta_mixed:
        mixed_pos = []
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

    if n_counts_frame.shape[0] != 0:
        stripped_frame, grouped_min_frame = get_grouping(n_counts_frame)
        final_frame = mixed_counts.merge(positions_frame, on='sample_name', how='left'). \
            merge(stripped_frame, on='sample_name', how='left').drop(['standard_name'], axis=1)
        final_frame['Exclude_from_analysis'] = np.where(final_frame['sample_name'].isin
                                                        (grouped_min_frame['sample_name']), "N", "Y")
        if args.metadata is not None:
            with_metadata = pd.merge(pd.read_csv(args.metadata).drop(['Exclude_from_analysis'], axis = 1), final_frame,
                                     how='inner', left_on='WGS_Id', right_on='sample_name').drop(['sample_name'],
                                                                                                 axis=1)
            with_metadata.sort_values(by=['WGS_Id']).to_csv("identical_sequences_w_metadata.csv", index=False)
        if args.metadata is None:
            final_frame.sort_values(by=['sample_name']).to_csv("identical_sequences.csv", index=False)


if __name__ == '__main__':
    main()

