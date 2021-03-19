# get the list of unique sample ids from the sam file directory
# glob wildcards will find the portion of all the filename in the directory based on the pattern below and save them to a tuple

# will ignore all other files of different extension
# NOTE: for glob_wildscards, ALWAYS need to put a comma after the variable for proper separation

ids, = glob_wildcards("fastq/{id}.fastq.gz")

# specify different wildcard for separating read 1 and 2 for paired read alignment
unique_id, = glob_wildcards("fastq/{read_id}_R1_001.fastq.gz")

# identify the unique sample id separate from the lane and read id that will be used for ivar

clean_ids, = glob_wildcards("fastq/{clean_id}_L001_R1_001.fastq.gz")


# can specify the output directory for certain command line utilities
# NOTE: do not specify a different output directory for the outputs of the ivar commands or the workflow will error out


# declare the global outputs from the different rules

# using comma and declaring variable names, multiple outputs can be declared for different rules
# need to use expand in the rule all inputs but not in the patterns for the individual rules
# can use a single variable in expand to fill the inputs with the output directory of choice 

rule all: 
	input: 
		html=expand("fastqc/{id}_fastqc.html", id=ids),
		zip=expand("fastqc/{id}_fastqc.zip", id=ids),
		bams=expand("bam/{read_id}_aligned.sorted.bam", read_id=unique_id),
		merged_bam=expand("ivar/{clean_id}_aligned.sorted.merged.bam", clean_id = clean_ids),
		index=expand("ivar/{clean_id}_aligned.sorted.merged.bam.bai", clean_id = clean_ids),
		trim_bam=expand("ivar/{clean_id}.primertrimmed.bam", clean_id = clean_ids),
		sorted_trim_bam=expand("ivar/{clean_id}.merged.sorted.bam", clean_id = clean_ids),
		fa=expand("ivar/{clean_id}.consensus.fa", clean_id = clean_ids),
		txt=expand("ivar/{clean_id}.consensus.qual.txt", clean_id = clean_ids),
		vars=expand("ivar/{clean_id}.variants.tsv", clean_id = clean_ids)
			

# create fastqc report of the raw fastqs for QC records

rule fastqc: 
	input: 
		"fastq/{id}.fastq.gz"

	output: 
		"fastqc/{id}_fastqc.html",
		"fastqc/{id}_fastqc.zip"

	shell: 
		"fastqc {input} -o fastqc/"

# perform paired read alignment using bwa mem

rule align: 
	input: 
		read_1="fastq/{read_id}_R1_001.fastq.gz",
		read_2="fastq/{read_id}_R2_001.fastq.gz",
		ref="/home/mwatson/COVID-19/reference/MN908947.fa"

	output: 
		"bam/{read_id}_aligned.sorted.bam"
	shell: 
		"bwa mem -t 4 {input.ref} {input.read_1} {input.read_2} | samtools sort | samtools view -S -b > {output}"

# merge the unique samples from all lanes into one, theinput requires a specific line for each lane to merge based on unique sample id
# specify new output directory ivar for ivar outputs

rule merge: 
	input: 
		l1="bam/{clean_id}_L001_aligned.sorted.bam",
		l2="bam/{clean_id}_L002_aligned.sorted.bam",
		l3="bam/{clean_id}_L003_aligned.sorted.bam",
		l4="bam/{clean_id}_L004_aligned.sorted.bam"

	output: "ivar/{clean_id}_aligned.sorted.merged.bam"

	shell: 
		"samtools merge {output} {input.l1} {input.l2} {input.l3} {input.l4}"

# index the sorted bam file

rule index: 
	input: 
		"ivar/{clean_id}_aligned.sorted.merged.bam"

	output: 
		"ivar/{clean_id}_aligned.sorted.merged.bam.bai"

	shell: 
		"samtools index {input}"


# output the trimmed bam files to the same directory as sorted through prefix
# use params to include the wildcard in the shell command for the sample prefix requirement

# the prefix field for ivar can specify the output directory as well as the sample pattern name

rule trim: 
	input: 
		bam="ivar/{clean_id}_aligned.sorted.merged.bam",
		bed="/home/mwatson/COVID-19/reference/nCoV-2019.primer.bed"

	output: 
		"ivar/{clean_id}.primertrimmed.bam"	

	params: 
		sam="{clean_id}"

	shell: 
		"ivar trim -i {input.bam} -b {input.bed} -p ivar/{params.sam}.primertrimmed -q 15 -m 20"

# generate a consensus fasta from the trimmed bam file
# switch the consensus input to the non-trimmed bam for experimentation

rule sort_trim: 
	input: 
		"ivar/{clean_id}_aligned.sorted.merged.bam"

	output: 
		"ivar/{clean_id}.merged.sorted.bam"

	shell: 
		"samtools sort {input} > {output}"


rule consensus: 
	input: 
		"ivar/{clean_id}.merged.sorted.bam"
	output: 
		fa="ivar/{clean_id}.consensus.fa",
		txt="ivar/{clean_id}.consensus.qual.txt"
	params: 
		fa="{clean_id}"
	shell: 
		"samtools mpileup -d 0 -aa -A -Q 0 {input} | ivar consensus -p ivar/{params.fa}.consensus -q 20 -t 0.75"

rule variants: 
	input: 
		ref="/home/mwatson/COVID-19/reference/MN908947.fa",
		bam="ivar/{clean_id}.merged.sorted.bam"
	output: 
		"ivar/{clean_id}.variants.tsv"
	params: 
		id="{clean_id}"
	shell: 
		"samtools mpileup -aa -A -d 600000 -B -Q 0 {input.bam} | ivar variants -p ivar/{params.id}.variants -q 20 -t 0.03 -r {input.ref}"


		





	




