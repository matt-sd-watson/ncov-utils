import pandas as pd
import os
from pathlib import Path
import glob

# Author   Matthew Watson  May 4, 2021
# In progress


var_dir = "/NetDrive/Projects/COVID-19/BCC_analysis_results/Plates/"

var_dict = {}
processed = 0
# search in specific directories with variant extensions to find variant.tsv files
for x in glob.glob(var_dir + "*/*callVariants/"):
    for path in Path(x).rglob('*variants.tsv'):
        var_frame = pd.read_csv(os.path.abspath(path), sep='\t')
        vars_list = []
        # ensure the vars frame is not empty and that the WGS Id is uniquw
        if 'ALT_DP' in var_frame and len(var_frame.index) != 0 and os.path.basename(path).split(".")[0].split("_S")[0] \
                not in var_dict:
            processed += 1
            print(os.path.basename(path).split(".")[0].split("_S")[0] + " " + "(" + str(processed) + ")")
            for index, row in var_frame.iterrows():
                if row["ALT_DP"] >= 10:
                    vars_list.append(str(row["REF"]) + str(row["POS"]) + str(row["ALT"]))
            var_dict[os.path.basename(path).split(".")[0].split("_S")[0]] = ",".join(vars_list)

vars_data_frame = pd.DataFrame(var_dict.items(), columns=["sample", "variants"])
vars_data_frame.to_csv("all_vars.csv", index=False)
