from Bio import SeqIO
import argparse


def main():
    parser = argparse.ArgumentParser(description='Verify that a multi-FASTA has sufficiently long sequences for'
                                                 'Nextstrain phylodynamic analysis and proper genome completeness')
    parser.add_argument('--input_fasta', '-i', type=str, help='input multi-FASTA to verify', required=True)

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input_fasta), 'fasta')

    too_short = 0
    for record in fasta_sequences:
        if len(record) < 29000:
            too_short += 1
            print("WARNING: The following sequences are too short for analysis: Sequence: {}, Length: {}".format(record.name, len(record)))
        try:
            genome_completeness = 100 - (100 * record.seq.count("N") / len(record))
            if genome_completeness <= 90:
                print("WARNING: The following sequences have genome completeness below 90%: {}  {}".format(
                    record.name, genome_completeness))
        except ZeroDivisionError:
            print("WARNING: sequence does not have any length: {}".format(record.name) + "\n")
    if too_short == 0:
        print("All sequences are long enough for analysis")


if __name__ == '__main__':
    main()

