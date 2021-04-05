from os import listdir
from os.path import isfile, join
import argparse

"""
This script takes an input a directory containing output ncov-19 samples from the ivar pipeline
For PHO samples, the way to designate a repeat sample is to place -v# after the unique WGS_Id, where
# is the repeat number
The script is esigned to find all previous versions of the sample 
i.e. if sample-v3 is in the directory, the output will save sample, and sample-v2

output: 
	= txt file containing all the sample repeat names as well as the original and previous sample names. 
	The output txt file is suitable for use with faSomeRecords to subset the master multi-fasta
"""

parser = argparse.ArgumentParser(description='Identify the set of sample repeats from an input directory')
parser.add_argument('--input_dir', '-i', type=str, help='Directory containing ncov samples', required=True)

args = parser.parse_args()

files = [f for f in listdir(args.input_dir) if isfile(join(args.input_dir, f))]

repeat_names = set([i.split('_S')[0] for i in files if '-v' in i])

strip_number = [i.split('-v')[1] for i in repeat_names]
strip_dict = dict(zip(repeat_names, strip_number))

total_names = []

for key, value in strip_dict.items():
    total_names.append(key.split('-v')[0])
    if int(value) <= 2:
        total_names.append(str(key.split('-v')[0] + "-v" + str(2)))
    else:
        for i in range(2, int(value) + 1, 1):
            total_names.append(str(key.split('-v')[0] + "-v" + str(i)))

with open("repeats.txt", "w") as handle:
    for line in sorted(set(total_names)):
        handle.write(line + "\n")
