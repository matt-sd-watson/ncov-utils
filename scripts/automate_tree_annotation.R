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
if (args$lineage_id == "B.1.1.7") {
  random_remove <- sample_n(tr.df.labs, as.integer(0.85*nrow(tr.df.labs)))
} else {
  random_remove <- sample_n(tr.df.labs, 0)
}

tree <- drop.tip(tree, random_remove$label)

## Create a dataframe of the tree
tr.df <- fortify(tree)
## Create a list of tree labels
tr.df.labs <- tr.df %>%
  filter(isTip == "TRUE") %>%
  select(label)

meta.df <- tr.df.labs %>%
  mutate(label.2 = ifelse(grepl("PHLON2", label), "ON-PHO", "non-PHO"))

# treat the countries as a factor for plotting
meta.df$label.2 <- as.factor(meta.df$label.2)

## Define colours for tree annotation
#3 must have the same name annotations as the label.2 created above
cols <- c("ON-PHO" = "#08E8DE",  "non-PHO" = "dark blue")

tr.df$nml_lab <- paste("ON-PHL", str_split_fixed(tr.df$label, "PHLON|-SARS", 4)[,2],
                       str_split_fixed(tr.df$label, "PHLON|-SARS", 4)[,3], sep = "-")

meta.df <- merge(meta.df, tr.df, by.x = "label", by.y = "label", all = TRUE)

pl.1 <- ggtree(tree, size = 0.25)

pl.2 <- pl.1 %<+% meta.df +
  geom_tippoint(aes(x=x+0.000001, subset=label.2 == "ON-PHO", label = label, colour = label.2),  size = 1.2, shape = 16) +
  geom_tippoint(aes(x=x+0.000001, subset=label.2 == "non-PHO", label = label.2, colour = label.2),  size = 1.2, shape = 16) +  
  geom_tiplab(aes(x=x+0.000001, subset=label.2 == "ON-PHO", label = nml_lab),  size = 1.5, offset = 0.000005) +
  geom_tiplab(aes(x=x+0.000001, subset=label.2 == "non-PHO", label = label),  size = 1.5, offset = 0.000005) +
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
pdf(paste(paste(args$lineage_id, "_phylo_tree_", sep = ""), format(Sys.Date(), '%d%b%Y'), ".pdf", sep=""), 
    width = 11, height = 15)
pl.2
dev.off()





