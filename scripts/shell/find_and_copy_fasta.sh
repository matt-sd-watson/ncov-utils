# this script takes an input a list of filenames, a directory to search in, and a final output directory
# It will find all of the names in the input list from the first directory, and copy them into the output directory

# read the input file, find the fasta with the name, and copy to target directory
cat $1 | while read line 
do
  find $2 -type f \( -name "$line*.fa" -o -name "$line*.fasta" \) -exec cp {} $3 \;
done

# Print the number of input lines and output files copied to ensure match
input_lines=$(less $1 | wc -l)
output_size=$(ls $3*.fa | wc -l)

echo "Number of target files:$input_lines"
echo "Number of files copied:$output_size"






