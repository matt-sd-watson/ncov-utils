
start=`date +%s`

rm -rf splits/ temp_lin/

mkdir -p splits

cp $1 splits/

round() { awk -v n=$1 -v d=$2 'BEGIN{print int((n+d/2)/d) * d}'; }

# automatically split the complete fasta into chunks of 1500

let intervals=$(grep ">" $1 | wc -l)/1500

pyfasta split -n $(round $intervals 1) splits/$1

rm splits/$1

source ~/anaconda3/etc/profile.d/conda.sh
conda activate pangolin

mkdir -p temp_lin

ls splits/*.fa | xargs -n 1 basename | parallel -j $(nproc) 'pangolin splits/{} --outfile temp_lin/{}.csv'

awk '(NR == 1) || (FNR > 1)' temp_lin/*.csv > lineage_report_cat.csv

rm -r splits/ temp_lin/

end=`date +%s`
runtime=$((end-start))

echo "Running pangolin took $runtime seconds"



