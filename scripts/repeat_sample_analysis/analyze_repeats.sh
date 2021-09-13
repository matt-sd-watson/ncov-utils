
# positional argument #1 = .txt file with the names of all samples for inspection
# positional argument #1 - output directory
# optional positional argument #3 - metadata from the database to merge with the repeat analysis
# example command - sh analyse_repeats.sh names.txt . 

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P) || exit

mkdir -p $2

python $dir/faSomeRecords.py -f /NetDrive/Projects/COVID-19/Other/master_fasta/complete*.fa -l $1 -o $2/repeats.fa

augur align \
  --sequences $2/repeats.fa \
   --nthreads auto \
   --reference-sequence $dir/reference.gb \
   --output $2/repeats_aln.fasta


snp-sites -c -o $2/repeats_snps_only.fa $2/repeats_aln.fasta

snp-sites -o $2/repeats_snps_only_mixed.fa $2/repeats_aln.fasta

snp-dists -m -c $2/repeats_aln.fasta > $2/repeats_snp_dists.csv

if [ ! -z $3 ] 
then 
    python $dir/analyze_repeats.py -f $2/repeats.fa -s $2/repeats_snp_dists.csv -o $2/ -m $2/repeats_snps_only_mixed.fa -d $3
else
    python $dir/analyze_repeats.py -f $2/repeats.fa -s $2/repeats_snp_dists.csv -o $2/ -m $2/repeats_snps_only_mixed.fa
fi






