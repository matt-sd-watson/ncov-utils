
mkdir -p $2

/usr/local/bin/./faSomeRecords /NetDrive/Projects/COVID-19/Other/master_fasta/complete*.fa $1 $2/repeats.fa

augur align \
  --sequences $2/repeats.fa \
   --nthreads auto \
   --reference-sequence /home/mwatson/COVID-19/reference/reference.gb \
   --output $2/repeats_aln.fasta


snp-sites -c -o $2/repeats_snps_only.fa $2/repeats_aln.fasta

snp-sites -o $2/repeats_snps_only_mixed.fa $2/repeats_aln.fasta

snp-dists -m -c $2/repeats_aln.fasta > $2/repeats_snp_dists.csv

if [ ! -z $3 ] 
then 
    python /home/mwatson/COVID-19/utilities/analyze_repeats.py -f $2/repeats.fa -s $2/repeats_snp_dists.csv -o $2 -m $2/repeats_snps_only_mixed.fa -d $3
else
    python /home/mwatson/COVID-19/utilities/analyze_repeats.py -f $2/repeats.fa -s $2/repeats_snp_dists.csv -o $2 -m $2/repeats_snps_only_mixed.fa
fi




