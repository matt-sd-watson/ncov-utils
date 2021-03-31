library(tidyverse)
#install.packages("ggtree")
library(ggtree)
# install.packages("data.table")
library(data.table)
# install.packages("phytools")
library(phytools)
library(ape)
library(grid)
library(stringr)
library(Biostrings)
library(TraMineR)
library(ggplot2)
library(argparse)
library(lubridate)

## Inputs: 
## - a tree in Newick format generated through iqtree
## - A lineage ID corresponding to the Pango Lineage of the tree samples
## A qc50 list with the names of samples that pass a QC threshold as "WGS_Id"

## Outputs: 
## - A customized annotated tree in pdf format for each the lineage samples

parser <- ArgumentParser(description='Annotate a phylogenetic tree generated through iqtree with custom labelling')
parser$add_argument('--input_tree', type = "character", 
                    help='path to the tree input')
parser$add_argument('--output_directory', type = "character", 
                    help='output directory for annotated tree in a PDF report')
parser$add_argument('--lineage_id', type = "character", 
                    help='Pango lineage identifier for file output annotation')
parser$add_argument('--qc50_list', type = "character", 
                    help='Metdata for samples that pass QC50 threshold')

args <- parser$parse_args()

tree <- read.tree(args$input_tree)
qc50 <- read.table(args$qc50_list, header = T, sep = ',', fill = TRUE, quote = "")

tree <- drop.tip(tree, "MN908947")

## Create a dataframe of the tree
tr.df <- fortify(tree)
## Create a list of tree labels
tr.df.labs <- tr.df %>%
  filter(isTip == "TRUE") %>%
  select(label)

# if the sample is not on the qc50 list, drop it from the tree
sample_no_pass <- tr.df.labs[!tr.df.labs$label %in% qc50$WGS_Id,]
tree <- drop.tip(tree, sample_no_pass$label)

# sub sample for B.1.1.7 because the number of samples is too large to view in a tree

## Create a dataframe of the tree
tr.df <- fortify(tree)
## Create a list of tree labels
tr.df.labs <- tr.df %>%
  filter(isTip == "TRUE") %>%
  select(label)

qc50$Date <- as.Date(qc50$upload_date, format = "%m/%d/%Y")
date_cat <- interval(ymd(Sys.Date() - 7), ymd(Sys.Date()))

qc50$date_cat <- ifelse(qc50$Date %within% date_cat, "Latest Week", "All Past Weeks")

qc50_keep <- subset(qc50, select = c(WGS_Id, Date, date_cat))

meta.df <- merge(tr.df.labs, qc50_keep, by.x = "label", by.y = "WGS_Id")

# treat the countries as a factor for plotting
meta.df$date_cat <- as.factor(meta.df$date_cat)

# if the lineage is B.1.1.7, set the background tips to a transparent grey
# otherwise, have a darker grey for lineage trees with fewer samples
grey_trans <- rgb(179, 179, 179, max = 255, alpha = 45)

cols <- c()
## Define colours for tree annotation depending on lineage identifier
if (args$lineage_id == "B.1.1.7") {
  cols <- c("Latest Week" = "dark blue",  "All Past Weeks" = grey_trans)
} else {
  cols <- c("Latest Week" = "dark blue",  "All Past Weeks" = "grey60")
}

meta.df$nml_lab <- paste("ON-PHL", str_split_fixed(meta.df$label, "PHLON|-SARS", 4)[,2],
                         str_split_fixed(meta.df$label, "PHLON|-SARS", 4)[,3], sep = "-")

pl.1 <- ggtree(tree, size = 0.25)

pl.2 <- pl.1 %<+% meta.df +
  geom_tippoint(aes(x=x+0.000001, subset=date_cat == "Latest Week", label = label, colour = date_cat),  size = 1.2, shape = 16) +
  geom_tippoint(aes(x=x+0.000001, subset=date_cat == "All Past Weeks", label = label, colour = date_cat),  size = 1.2, shape = 16) +  
  # geom_tiplab(aes(x=x+0.000001, subset=label.2 == "ON-PHO", label = nml_lab),  size = 1.5, offset = 0.000005) +
  # geom_tiplab(aes(x=x+0.000001, subset=label.2 == "non-PHO", label = label),  size = 1.5, offset = 0.000005) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        legend.text=element_text(size=6),
        legend.position = c(0.85, 0.6),
        legend.title=element_text(size=7, face = "bold"),
        legend.box = "vertical",
        legend.box.margin = margin(0.01,0.02,0.02,0.02,"cm"),
        legend.box.background = element_rect(colour = "grey50"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.key.size = unit(3,"lines")) +
  scale_color_manual("Location", values=cols) +
  guides(color = guide_legend(override.aes = list(size = 1.75)))

setwd(args$output_directory)

pl.2
## SAVE annotated tree

plot_name <- paste(paste(args$lineage_id, "_phylo_tree_", sep = ""), format(Sys.Date(), '%d%b%Y'), ".pdf", sep="")

ggsave(filename = plot_name, plot = pl.2, width = 11, height = 15)

dev.off()
rm(list = ls())






