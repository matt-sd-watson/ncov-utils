# Author- Matthew Watson	January 18, 2021

pangolin $1 --outfile lineage_report.csv

# this script adds a new column to the pangolin
# lineage report with the pangolin version number when it was run (i.e. v2.1.6)
# to facilitate better exporting into the Access database set up by Jen Guthrie


# make a variable using the pangolin version to input into csv
version=$(pangolin --version | sed 's/^pangolin //')
final_version="v$version"


# make the file tab sep for later header formatting
sed 's/,/\t/g' lineage_report.csv > lineage_report.tsv

# make a new column with the version variable
# this assumes that the incoming pangolin report has only 6 columns so insert at 7th

# insert pangolin version into the last column and re-add the headers

awk -v var="$final_version" 'BEGIN{FS=OFS="\t"}{$7=$7" "var; print }' lineage_report.tsv | awk 'BEGIN {FS=OFS="\t"; OFS=","} {print $1, $2, $3, $4, $5, $7, $6}' | sed -n '1!p' | sed -e '1i\taxon\tlineage\tprobability\tpangoLEARN_version\tstatus\tpangolin_version\tnote' > $2

sed -i 's/\t/,/g' $2


# maintain the tab separation for the output file
# perl -p -i -e 's/ /\t/g' output.tsv

# remove the first header and replace with header that includes the proper column names
# switch the 6th and 7th columns to have the unstructured notes portion as the rightmost column for formatting purposes
# switch around the columns

# sed -n '1!p'  lineage_report_2.tsv | sed -e '1i\taxon\tlineage\tprobability\tpangoLEARN_version\tstatus\tnote\tpangolin_version' > $2

# cut -f1-9 output_2.csv > $2

rm lineage_report.csv
#rm output.tsv
rm lineage_report.tsv
#rm output_2.csv




