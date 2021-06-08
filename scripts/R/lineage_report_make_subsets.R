library(argparse)
library(dplyr)
library(data.table)

parser <- ArgumentParser(description='Create subsets for lineage reporting')
parser$add_argument('--input_gisaid', type = "character", 
                    help='Input metadata CSV for gisaid context sequences')
parser$add_argument('--input_lineage', type = "character", 
                    help='Input lineage report CSV for PHO sequences')
parser$add_argument('--output_filenames', type = "character", 
                    help='list of PHO sample names to subset from')
parser$add_argument('--output_gisaid_names', type = "character", 
                    help='list of Gisaid context sample names to subset from')
parser$add_argument('--lineage_id', type = "character", 
                    help='ID of the lineage to create subset')
parser$add_argument('--subset_number', type = "integer", 
                    help='Number of Sequences to randomly subset from Gisaid')
parser$add_argument('--metadata_output', type = "character", 
                    help='Metadata output for Gisaid subset sequences')

# parser$add_argument('--output_dir', type = "character", 
                    # help='target output dir')

args <- parser$parse_args()

# read in a metadata table from gisaid and filter based on criteria for the lineage report
metadata <- read.table(args$input_gisaid, header = TRUE, sep='\t', fill = TRUE, quote = "")
dim(metadata)

# format Gisaid date and remove any duplicates
metadata$date <- as.Date(metadata$date, format = "%Y-%m-%d")

metadata_complete <- metadata[complete.cases(metadata$date),]

# make the subset without Canadian samples and add them later
subset_gisaid <- subset(metadata_complete, pango_lineage == args$lineage_id & country != "Canada")

# add in all Canadian samples after sub-sampling
subset_to_add <- subset(metadata_complete, pango_lineage == args$lineage_id & country == "Canada")
# exclude the Gisaid samples that came from PHO (have ON-PHL in the strain name)
subset_to_add <- subset_to_add[!grepl("ON-PHL", subset_to_add$strain),]

print("Number of non-PHO Canadian sequences:", quote = F)
print(nrow(subset_to_add))

#calculate the difference between the requested subset number
# and the number of possible records to include
difference <- nrow(subset_to_add) + nrow(subset_gisaid) - as.numeric(args$subset_number)

# if the requested number is more than total that are available, keep all sequences
# otherwise, determine how to top up
if (difference <= 0) {
  subset_gisaid <- subset_gisaid
  subset_to_add <- subset_to_add
} else {
  difference_top_up <- as.numeric(args$subset_number) - nrow(subset_to_add)
  # if the number of canadian sequences is greater than the requested amount, do not
  # add any more sequences
  # otherwise, top up the total with non canadian sequences
  if (difference_top_up <= 0) {
    subset_gisaid <- sample_n(subset_gisaid, 0)
    subset_to_add <- subset_to_add
  }  else {
    subset_gisaid <- sample_n(subset_gisaid, difference_top_up)
    subset_to_add <- subset_to_add
  }
}

# create final subset with random and Canadian sequences
final_subset <- rbind(subset_gisaid, subset_to_add)

setnames(final_subset, old = c("strain", "date", "gisaid_epi_isl", "Nextstrain_clade",
                                  "pango_lineage", "GISAID_clade"),
         new = c("sequence_name", "sample_date", "central_sample_id","phylotype",
                 "lineage", "uk_lineage"))

# setwd(args$output_dir)
# write the samples names to a txt file to sub-sample with faSomeRecords
write.table((unique(final_subset$sequence_name)), file = args$output_gisaid_names, sep="\t", col.names = F, row.names = F,
            quote = F)

# establish the final column for writing the metadata to avoid column spill over
end_columns <- as.numeric(which( colnames(final_subset)=="uk_lineage" ))

write.csv(final_subset[,1:end_columns], file = args$metadata_output, row.names = F,
            quote = F)

# read the pangolin lineage report and get the WGS Ids with the lineage
pangolin_pho <- read.table(args$input_lineage, header = TRUE, sep=',', fill = TRUE, quote = "")

subset_pho <- subset(pangolin_pho, lineage == args$lineage_id)

# write the unique IDs only
unique_ids <- unique(subset_pho$taxon)

write.table(unique_ids, file = args$output_filenames, sep="\t", col.names = F, row.names = F,
            quote = F)