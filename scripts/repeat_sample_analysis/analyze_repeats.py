import pandas as pd
from Bio import SeqIO
import argparse
import numpy as np
import warnings

warnings.filterwarnings("ignore")


# return a frame that selects the sample with the highest completeness among repeats
def get_lowest_n_count_grouping(master_frame):
    # master_frame['standard_name'] = master_frame['sample_name'].map(lambda x: x.split('-v', 1)[0])
    min_frame = master_frame.loc[master_frame.groupby('standard_name').N_counts.idxmin()]

    return min_frame


# create a data frame thatr groups samples based on their shared name, creating repeat groups
def get_grouped_data_frame(data_frame):
    data_frame['standard_name'] = data_frame['sample_name'].map(lambda x: x.split('-v', 1)[0])
    return data_frame


# filter a data frame to exclude samples that are contained in a second data frame
def filter_first_frame(data_frame_1, data_frame_2):
    return data_frame_1[(~data_frame_1.standard_name.isin(data_frame_2.sample_name)) &
                        (~data_frame_1.sample_name.isin(data_frame_2.standard_name))]


# create a dataframe that contains only repeat groups where all samples are above 90% completeness
# select the one with the lowest number of mixed counts. in the case of a tie, select the one with higher
# completeness
def whole_group_over_90(data_frame):
    all_over_90 = data_frame.groupby('standard_name').filter(
        lambda x: (x['N_counts'] <= 2990).all())
    return all_over_90.loc[
        all_over_90.sort_values(by=['N_counts'], ascending=True).groupby('standard_name').mixed_counts.idxmin()]


def main():
    parser = argparse.ArgumentParser(description='Analyze a set of sample repeats for similarity and'
                                                 'coverage statistics')
    parser.add_argument('--multi_fasta', '-f', type=str, help='Multi-fasta to evaluate sequences from', required=True)
    parser.add_argument('--snp_dists', '-s', type=str, help='CSV file of SNP differences from snp-dists in'
                                                            'molten format', required=True)
    parser.add_argument('--output_dir', '-o', type=str, help='output directory for single modified fasta files',
                        required=True)
    parser.add_argument('--fasta_mixed', '-m', type=str, help='Multi-fasta without exclusively ACGT to count'
                                                              'mixed positions', required=True)
    parser.add_argument('--metadata', '-d', type=str, help='Optional metadata frame with additional sample information',
                        required=False)

    args = parser.parse_args()

    snp_dists = pd.read_csv(args.snp_dists, header=None)
    snp_dists.columns = ['sam_1', 'comparing', 'distance']

    # filter if the sample 1 name is in the second column, and filter by snp distances
    new_frame = snp_dists[[x[0] in x[1] for x in zip(snp_dists['sam_1'], snp_dists['comparing'])]]
    identical_frame = new_frame.loc[(new_frame['sam_1'] != new_frame['comparing'])][snp_dists.distance == 0]
    non_identical_frame = new_frame.loc[(new_frame['sam_1'] != new_frame['comparing'])][snp_dists.distance != 0]

    diff_name_identical = snp_dists[~snp_dists.sam_1.isin(identical_frame.sam_1)]
    diff_name_identical = diff_name_identical[~diff_name_identical.sam_1.isin(identical_frame.comparing)]
    diff_name_identical = diff_name_identical[snp_dists.apply(lambda x: x.sam_1 not in x.comparing,
                                                              axis=1)][snp_dists.distance == 0]
    diff_name_identical = diff_name_identical[~diff_name_identical['sam_1'].isin(['MN908947'])]

    if len(diff_name_identical) != 0:
        diff_name_identical.to_csv("different_names_identical.csv", index=False)

    fasta_sequences = SeqIO.parse(open(args.multi_fasta), 'fasta')

    sam_n_counts = []
    non_identical = []
    for record in fasta_sequences:
        keys = ['sample_name', 'N_counts']
        values = [record.id, record.seq.count("N")]
        to_add = dict(zip(keys, values))
        sam_n_counts.append(to_add)
        if record.id in non_identical_frame.sam_1.unique() or record.id in non_identical_frame.comparing.unique():
            non_identical.append(record.id)

    n_counts_frame = pd.DataFrame(sam_n_counts)
    non_identical.sort()

    print("Identical Sequences: {}".format(n_counts_frame.shape[0]))
    if n_counts_frame.shape[0] != 0:
        print(n_counts_frame.sort_values(by=['sample_name'], ascending=True))
    else:
        print("None")

    print("Non-Identical sequences: {}".format(len(non_identical)))
    if len(non_identical) != 0:
        with open("Non-identical sequences.txt", "w") as handle:
            for lines in non_identical:
                handle.write(lines+"\n")
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
        count = 0
        for index, item in enumerate(record.seq):
            if item not in nuc:
                count += 1
                mixed_pos.append("{}:{}".format(item, int(index + 1)))

        keys = ['sample_name', 'mixed_counts']
        if record.id != "MN908947":
            values = [str(record.id), count]
            if values:
                to_add = dict(zip(keys, values))
                mixed_post.append(to_add)
                pos[str(record.id)] = ';'.join(n.strip("[]") for n in mixed_pos)
            else:
                raise TypeError("The values are undefined")

    mixed_counts = pd.DataFrame(mixed_post)
    mixed_counts["identical_called_snps"] = np.where(mixed_counts['sample_name'].isin(non_identical), "N", "Y")
    positions_frame = pd.DataFrame(pos.items(), columns=['sample_name', 'mixed_positions'])

    if n_counts_frame.shape[0] != 0:
        final_frame = mixed_counts.merge(positions_frame, on='sample_name', how='left').merge(n_counts_frame,
                                                                                              on="sample_name",
                                                                                              how='left')

        # filter all repeat groups where the coverage is over 90%
        final_frame_over_90 = whole_group_over_90(get_grouped_data_frame(final_frame))

        # get thw lowest N counts for samples that were not in the both 90% category
        # if a sample of any of its partners is in the over 90% category, do not include them in the
        # N count lowest frame
        min_frame_unfiltered = get_lowest_n_count_grouping(get_grouped_data_frame(final_frame))
        min_frame_filtered = filter_first_frame(min_frame_unfiltered, final_frame_over_90)

        # if the sample is in either of the minimum frames, do not exclude
        final_frame["Exclude_from_analysis"] = np.where((final_frame["sample_name"].isin
                                                         (min_frame_filtered["sample_name"])) |
                                                        (final_frame["sample_name"].isin
                                                         (final_frame_over_90["sample_name"])), "N", "Y")

        # scenarios:
        # 1. if SNPs match, it has fewer than 2 mixed SNPs, and it is excluded from analysis, teg:
        #       "repeat- no SNP mismatch"
        # 2. if SNPs match, it has more than 2 mixed SNPs, and it is excluded from analysis, teg:
        #       "repeat- SNP mismatch- more than 2 mixed sites"
        # 3. if SNPS do not match, tag:
        #       "repeat- SNP mismatch - sample mix-up"
        # 4. Otherwise, tag "Do not exclude"
        final_frame['Exclude_from_analysis_details'] = np.where((final_frame['identical_called_snps'] == "Y") &
                                                                (final_frame['mixed_counts'] <= 2) &
                                                                (final_frame["Exclude_from_analysis"] == "Y"),
                                                                "repeat- no SNP mismatch",
                                                                np.where((final_frame['identical_called_snps'] == "Y") &
                                                                         (final_frame['mixed_counts'] > 2) &
                                                                         (final_frame["Exclude_from_analysis"] == "Y"),
                                                                         "repeat- SNP mismatch- more than 2 "
                                                                         "mixed sites different",
                                                                         np.where((final_frame['identical_called_snps']
                                                                                   == "N") &
                                                                                  (final_frame['mixed_counts'] <= 2),
                                                                                  "repeat- SNP mismatch - "
                                                                                  "non identical",
                                                                                  np.where((final_frame
                                                                                            ['identical_called_snps'] ==
                                                                                            "N") &
                                                                                           (final_frame['mixed_counts']
                                                                                            > 2),
                                                                                           "repeat- SNP mismatch -"
                                                                                           "sample mix-up", ""))))

        if args.metadata is not None:
            with_metadata = pd.merge(pd.read_csv(args.metadata), final_frame, how='inner', left_on='WGS_Id',
                                     right_on='sample_name').drop(['sample_name'], axis=1)
            with_metadata.sort_values(by=['WGS_Id']).to_csv("identical_sequences_w_metadata.csv", index=False)
        else:
            final_frame.sort_values(by=['sample_name']).drop(['standard_name'], axis=1).to_csv("all_repeats.csv", index=False)


if __name__ == '__main__':
    main()

