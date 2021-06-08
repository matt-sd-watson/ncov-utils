import pysam
import collections
import re
import csv
import argparse

# Author    Matthew Watson  February 2, 2021
# This script accepts as input an indexed bam file of ncov19 from an ivar pipeline (Illumina)
# and a csv file of mutations in the first column in standard type variants format
# Note that only nucleotide-formatted mutations are accepted currently (i.e. no aa notations)
# It will parse the bam file and save the frequency of the nucleotide pileups for each mutation
# For deletions, the script will evaluate the pileup at the nucleotide in the middle of the deletion


def main():
    parser = argparse.ArgumentParser(description='Parse an indexed bam file for mutation frequencies at set '
                                                 'coordinates')
    parser.add_argument('--mutation_ids', '-m', type=str, help='CSV file with mutations of interest in type '
                                                               'variants format', required=True)
    parser.add_argument('--input_bam', '-i', type=str, help='Input ncov bam file with index file in the same '
                                                            'directory', required=True)
    parser.add_argument('--output_file', '-o', type=str, help='Output csv with mutation frequencies', required=True)

    args = parser.parse_args()

    sam_file = pysam.Samfile(args.input_bam, "rb")

    mutations = []
    with open(args.mutation_ids, "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for lines in csv_reader:
            mutations.append(lines[0])

    positions = {}
# create dict of the mutation name from the csv and the int position for the pileup
    for m in mutations:
        if "del" in m:
            match = re.match(r"([a-z|:]+)([0-9|:]+)", m, re.I)
            # if there is a deletion, evaluate the position in the middle of the deletion (rounded)
            del_length = int(match.groups()[1].split(":")[1])
            spot_del = int(match.groups()[1].split(":")[0]) + int(del_length / 2)
            # use the middle of the deletion as the key for the mutation
            positions[spot_del] = m
        else:
            match = re.match(r"([a-z]+)([0-9]+)", m, re.I)
        if match:
            # if match, use the nucleotide position as a key and the mutation name as value
            positions[int(match.groups()[1])] = m

# create a dict for the mutation name and th frequencies of nucleotides or deletions at the position
    sample_dict = {}
    for i in positions.keys():
        base_dict = []
        # ignore overlaps allows all paired end reads overlapping a position to be counted
        # ignoring overlaps will set the pileup numbers to be the same as they are viewed in tablet
        # https://github.com/pysam-developers/pysam/issues/703
        for pileupcolumn in sam_file.pileup(None, i - 1, i, ignore_overlaps=False):
            for pileupread in pileupcolumn.pileups:
                # pysam uses 0-based coordinates while the bam file is written with 1-based coordinates
                if pileupcolumn.pos == i - 1:
                    if pileupread.query_position is not None:
                        # if there is a gap or deletion, the query position will be None
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_dict.append(base)
                    else:
                        base_dict.append("Del")

    # use the mutation name as the key and the dictionary of frequencies as the value (nested dict)
        sample_dict[positions.get(i)] = str(dict(collections.Counter(base_dict)))

# remove curly brackets from the value that is a nested dictionary
    cleaned_dict = {key.strip(): item.strip('{|}').replace("'", "") for key, item in sample_dict.items()}

    with open(args.output_file, 'w') as csv_file:
        headers = ["sample", "mutation", "frequencies"]
        writer = csv.writer(csv_file)
        writer.writerow(headers)

        for key, value in cleaned_dict.items():
            writer.writerow([str(args.input_bam), key, value])


if __name__ == '__main__':
    main()

