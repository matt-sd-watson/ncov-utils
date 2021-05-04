import pandas as pd
import os
from pathlib import Path
import glob
import argparse

# Author   Matthew Watson  May 4, 2021
# In progress
# examples execution to search multiple NetDrive directories: 
# python parse_variants_cli.py -i "/NetDrive/Projects/COVID-19/BCC_analysis_results/Plates/*/*callVariants/" -i "/NetDrive/Projects/COVID-19/Analysis_Results/ncov2019ArticNf_1.6_results/" -o test_vars.csv
# search paths must be wrapped in quotes to be recognized by glob


def main():
    parser = argparse.ArgumentParser(description='parse a directory for variants TSV files from the ivar pipeline')
    parser.add_argument('--input_pattern', '-i', type=str, help='Input pattern containing the directory and'
                                                                'glob-appropriate pattern to search', required=True,
                        nargs='*')
    parser.add_argument('--output_file', '-o', type=str, help='Output CSV file', required=True)

    args = parser.parse_args()
    var_dict = {}
    processed = 0

    # search in specific directories with variant extensions to find variant.tsv files
    for i in args.input_pattern:
        for x in glob.glob(str(i)):
            for path in Path(x).rglob('*variants.tsv'):
                var_frame = pd.read_csv(os.path.abspath(path), sep='\t')
                # ensure the vars frame is not empty and that the WGS Id is unique
                if 'ALT_DP' in var_frame and len(var_frame.index) != 0:
                    vars_list = []
                    processed += 1
                    print("Processing: " + os.path.basename(path).split(".")[0].split("_S")[0] + " " + "(" +
                          str(processed) + ")")
                    for index, row in var_frame.iterrows():
                        if row["ALT_DP"] >= 10:
                            if "-" in row["ALT"]:
                                vars_list.append("del:" + str(row["REF"]) + str(row["POS"]) + str(row["ALT"]))
                            elif "+" in row["ALT"]:
                                vars_list.append("ins:" + str(row["REF"]) + str(row["POS"]) + str(row["ALT"]))
                            else:
                                vars_list.append(str(row["REF"]) + str(row["POS"]) + str(row["ALT"]))
                    var_dict[os.path.basename(path).split(".")[0].split("_S")[0]] = ",".join(vars_list)

        vars_data_frame = pd.DataFrame(var_dict.items(), columns=["sample", "variants"])
        vars_data_frame.replace("", float("NaN"), inplace=True)
        vars_data_frame.dropna(subset=["variants"], inplace=True)
        vars_data_frame.to_csv(args.output_file, index=False)


if __name__ == '__main__':
    main()
