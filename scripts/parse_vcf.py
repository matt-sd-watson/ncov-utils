import io
import os
import pandas as pd
import vcf
import argparse

# use a function to parse the main columns


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'region'})


def main():

    parser = argparse.ArgumentParser(description='read a series of VCF files and create a variants CSV')
    parser.add_argument('--input_file', '-i', type=str, nargs='+',
                        help='Series of input vcf files to parse and consolidate', required=True)
    parser.add_argument('--output_csv', '-o', type=str, help='output CSV file containing VCF information',
                        default='consolidated.CSV')

    args = parser.parse_args()

    master_frame = pd.DataFrame()

    for i in args.input_file:
        vars_frame = read_vcf(i)
        # remove consolidated columns
        vars_frame = vars_frame.drop(['INFO', 'FORMAT', 'unknown'], axis=1)

        # add column in the first position for sample name
        vars_frame.insert(0, column="WGS_Id", value=os.path.basename(i).split("_S")[0])
        # establish a list of extra columns to add
        columns_to_add = ['DP', 'DPB', 'RO', 'AO', 'PRO', 'PAO']

        # determine the column position where the new columns will be added
        column_index = len(vars_frame.columns)
        for value in columns_to_add:
            # increment the column position to continuously add desired columns
            column_index += 1
            vcf_reader = vcf.Reader(open(i, 'r'))
            values = []
            for record in vcf_reader:
                values.append(str(record.INFO[value]).strip("[|]"))
            vars_frame.insert(column_index - 1, column=value, value=values)
        # continue to merge the data frames for each vcf passed
        master_frame = pd.concat([master_frame, vars_frame])

    master_frame.to_csv(args.output_csv, index=False)


if __name__ == '__main__':
    main()
