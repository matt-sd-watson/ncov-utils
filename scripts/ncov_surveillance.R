library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(lubridate)
library(viridis)
library(gridExtra)
library(wesanderson)
library(RColorBrewer)
library(directlabels)
library(argparse)
library(ggpubr)

parser <- ArgumentParser(description='Process surveillance plots for bulk samples')
parser$add_argument('--input_qc', type = "character", 
                    help='Input QC CSV')
parser$add_argument('--input_nextclade', type = "character", 
                    help='Input Nextclade TSV')
parser$add_argument('--lineage_watchlist', type = "character", 
                    help='Input txt file of lineages to monitor')
parser$add_argument('--mutation_watchlist', type = "character", 
                    help='Input txt file of mutations to monitor')

# parser$add_argument('--output_dir', type = "character", 
# help='target output dir')

args <- parser$parse_args()

qc_data <- read.table(args$input_qc, header = T, sep = ',',
                      fill = TRUE, quote = "")

qc_data$date <- as.Date(qc_data$upload_date, format = "%m/%d/%Y")

# arrange by ascending date and remove all duplicate sample names to keep earliest database instance
qc_data <- qc_data %>% arrange(date)

qc_data <- qc_data[!duplicated(qc_data$WGS_Id),]

# read in the nextclade data that has the mutation profiles
nextclade_data <- read.table(args$input_nextclade, header = T, sep = '\t',
                             fill = TRUE, quote = "")

lineages <- as.vector(scan(args$lineage_watchlist, character(), quote = ""))

subset_lineages <- subset(qc_data, PANGO_lineage_updated %in% lineages) %>%
  arrange(date) %>%
  group_by(PANGO_lineage_updated, date) %>% 
  summarise(counts = n())

rolling_cases <- subset_lineages %>% group_by(PANGO_lineage_updated, .drop = FALSE) %>%
  mutate(cumcases = cumsum(counts))

setnames(rolling_cases, old = c('PANGO_lineage_updated'), new = c('lineage'))

plot_1 <- ggplot(rolling_cases, aes(x = date, y = cumcases)) + geom_point(col = "blue") + 
  facet_wrap(~lineage, ncol=2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0)) +
  geom_line(color='red',data = rolling_cases, aes(x = date, y = cumcases)) +
  geom_label_repel(aes(label = cumcases),
                   box.padding   = 0.2, 
                   point.padding = 0.3,
                   segment.color = 'grey50',
                   size = 3) +
  theme(axis.title.x = element_text(vjust=-0.5),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Variants of Concern (VOC)") + xlab("Logging Date") +
  ylab("Cumulative Case Counts") + scale_x_date(date_labels = "%b-%Y")

ggsave("cumulative_lineage_counts", plot = plot_1, device = "pdf", width=11, height=8.5)

mutations <- as.vector(scan(args$mutation_watchlist, character(), quote = ""))

# make sure that the nextclade data is taken from the qc50 sequences and not all
nextclade_qc50 <- nextclade_data[nextclade_data$seqName %in% qc_data$WGS_Id,]

# search for all mutations using grepl
nextclade_with_mut <- nextclade_qc50[grepl(paste(mutations,collapse="|"), 
                                           nextclade_qc50$aaSubstitutions),]

# remove any object matching pattern for concatenation
rm(list=ls(pattern="^mutation_frame.*"))

# create a mutation frame for each mutation with the samples containing the mutation
for (i in mutations) {
  nextclade_with_mut <- nextclade_qc50[grepl(i, nextclade_qc50$aaSubstitutions),]
  nextclade_with_mut$lineage <- i
  name <- paste("mutation_frame_", i, sep="")
  assign(name, nextclade_with_mut)
}

# bind all mutation frames
df_list <- mget(ls(pattern = "^mutation_frame.*"))
all_muts <- plyr::rbind.fill(df_list)

muts_frame_keep <- subset(all_muts, select = c(seqName, lineage))

mut_merged <- merge(qc_data, muts_frame_keep, by.x = "WGS_Id",
                    by.y = "seqName")

muts_rolling <- mut_merged %>% arrange(date) %>%
  group_by(lineage, date) %>%
  summarise(counts = n())

rolling_cases_muts <- muts_rolling %>% group_by(lineage) %>%
  mutate(cumcases = cumsum(counts))

plot_2 <- ggplot(rolling_cases_muts, aes(x = date, y = cumcases)) + geom_point(col = "blue") + 
  facet_wrap(~lineage) +
  theme_bw() +
  geom_line(color='red',data = rolling_cases_muts, aes(x = date, y = cumcases)) +
  geom_label_repel(aes(label = cumcases),
                   box.padding   = 0.2, 
                   point.padding = 0.3,
                   segment.color = 'grey50',
                   size = 3)  +
  theme(axis.title.x = element_text(vjust=-0.5),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Mutations of Interest") + xlab("Logging Date") +
  ylab("Cumulative Case Counts") + scale_x_date(date_labels = "%b-%Y")

ggsave("cumulative_mutation_counts", plot = plot_2, device = "pdf", width=11, height=8.5)

rolling_cases$designation <- "lineage"
rolling_cases_muts$designation <- "mutation"

date_cat <- interval(ymd(Sys.Date() - 7), ymd(Sys.Date()))

all_to_monitor <- rbind(rolling_cases, rolling_cases_muts)
all_to_monitor$time_cat <- ifelse(all_to_monitor$date %within% date_cat, "Latest Week", "All Past Weeks")

grouped_lin <- all_to_monitor[all_to_monitor$designation == "lineage",] %>%
  group_by(lineage, time_cat) %>% summarise(cat_counts = sum(counts))

plot_3 <- ggplot(grouped_lin, aes(x = time_cat, y = cat_counts, fill = time_cat)) + geom_bar(stat = 'identity') +
  facet_wrap(~lineage) + geom_text(aes(label=cat_counts), position=position_dodge(width=0.9), vjust=-0.25) +
  xlab("Time Category") + ylab("Case Counts") + scale_fill_manual("Time Category",
                                                                  values = c("dark blue", "dark red")) +
  ggtitle(paste("Variants of Concern (VOC) as of ", format(Sys.Date(), "%b-%d-%Y"), sep="")) +
  ylim(c(0, max(grouped_lin$cat_counts) + 35)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggsave("VOC_lineages_by_week", plot = plot_3, device = "pdf", width=11, height=8.5)

grouped_mut <- all_to_monitor[all_to_monitor$designation == "mutation",] %>%
  group_by(lineage, time_cat) %>% summarise(cat_counts = sum(counts))

plot_4 <- ggplot(grouped_mut, aes(x = time_cat, y = cat_counts, fill = time_cat)) + geom_bar(stat = 'identity') +
  facet_wrap(~lineage) + geom_text(aes(label=cat_counts), position=position_dodge(width=0.9), vjust=-0.25) +
  xlab("Time Category") + ylab("Case Counts") + scale_fill_manual("Time Category",
                                                                  values = c("dark blue", "dark red")) +
  ggtitle(paste("Mutations of Interest as of ", format(Sys.Date(), "%b-%d-%Y")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black")))

ggsave("Mutations_of_interest_by_week", plot = plot_4, device = "pdf", width=11, height=8.5)


colnames(muts_frame_keep) <- c("sample", "mutation")

all_lin_and_mut <- merge(muts_frame_keep, qc_data, by.x = "sample",
                         by.y = "WGS_Id") %>% select(sample, PANGO_lineage_updated,
                                                     mutation) %>% group_by(mutation,
                                                                            PANGO_lineage_updated) %>%
  summarise(counts = n(), percent = counts/sum(counts))

all_lin_and_mut$PANGO_lineage_updated <- as.factor(all_lin_and_mut$PANGO_lineage_updated)

mut_by_line_1 <- ggplot(all_lin_and_mut, aes(x = mutation, y = counts, fill = PANGO_lineage_updated)) +
  geom_bar(position="dodge", stat="identity") +
  xlab("AA Mutation") + ylab("Case Counts") + labs(fill='PANGOLIN lineage') +
  geom_text(aes(label=counts), position = position_dodge(width = 0.9), vjust = -0.25, size = 4) +
  ggtitle("Mutations of Interest") + ylim(c(0, max(all_lin_and_mut$counts) + 50)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

mut_by_line_2 <- ggplot(all_lin_and_mut, aes(x = mutation, y = counts, fill = PANGO_lineage_updated)) +
  geom_bar(position="fill", stat="identity") +
  xlab("AA Mutation") + ylab("Percentage of Case Counts") + labs(fill='PANGOLIN lineage') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_5 <- ggarrange(mut_by_line_1, mut_by_line_2, nrow = 2)

ggsave("Distribution_lineages_by_mutation", plot = plot_5, device = "pdf", width=11, height=8.5)

mut_lin_rolling_date <- subset(mut_merged, select = c(WGS_Id, PANGO_lineage_updated,
                                                      date, lineage))

colnames(mut_lin_rolling_date) <- c("name", "lineage", "date", "mutation")

mut_lin_date_grouped <- mut_lin_rolling_date %>% arrange(date) %>%
  group_by(mutation, lineage, date) %>%
  summarise(counts = n())

mut_lin_date_grouped_final <- mut_lin_date_grouped %>% group_by(mutation, lineage) %>%
  mutate(cumcases = cumsum(counts))

mut_lin_date_grouped_final_no_b117 <- mut_lin_date_grouped_final[mut_lin_date_grouped_final$lineage != "B.1.1.7",]

# get the most recent date for each of the lineages by mutation
max_date_each_lin_no_b117 <- mut_lin_date_grouped_final_no_b117 %>% group_by(lineage) %>% filter(date == max(date))

# plot the lineages separately for each mutation
# include a label at the most recent date of logging for each mutation
plot_6 <- ggplot(mut_lin_date_grouped_final_no_b117, aes(x = date, y = cumcases, col = lineage)) +
  geom_point(data = max_date_each_lin_no_b117, aes(color = lineage)) +
  geom_line() + facet_wrap(~mutation, ncol = 2) + ylab("Case Counts") +
  xlab("Logging Date") +
  theme(axis.title.x = element_text(vjust=-0.5),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(data = max_date_each_lin_no_b117,
                  aes(label = paste(lineage, paste("(", cumcases, ")", sep = ""), sep = " ")), hjust = 0, nudge_x = 0.05, size = 3,
                  col = "black") + ggtitle("Cumulative Counts of Lineages with Mutations of Interest, no B.1.1.7") +
  scale_x_date(date_labels = "%b-%Y")

mut_lin_date_grouped_final_b117 <- mut_lin_date_grouped_final[mut_lin_date_grouped_final$lineage == "B.1.1.7",]

# get the most recent date for each of the lineages by mutation to assign labels
max_date_each_lin_b117 <- mut_lin_date_grouped_final_b117 %>% group_by(lineage) %>% filter(date == max(date))

plot_7 <- ggplot(mut_lin_date_grouped_final_b117, aes(x = date, y = cumcases, col = lineage)) +
  geom_point(data = max_date_each_lin_b117, aes(color = lineage)) +
  geom_line() + facet_wrap(~mutation, ncol = 2) + ylab("Case Counts") +
  xlab("Logging Date") +
  theme(axis.title.x = element_text(vjust=-0.5),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(data = max_date_each_lin_b117,
                  aes(label = paste(lineage, paste("(", cumcases, ")", sep = ""), sep = " ")), hjust = 0, nudge_x = 0.05, size = 3,
                  col = "black") + ggtitle("Cumulative Counts of Lineages with Mutations of Interest, B.1.1.7") +
  scale_x_date(date_labels = "%b-%Y")

ggsave("cumulative_lineage_counts_watchlist_mutations_no_b117", plot = plot_6, device = "pdf", width=11, height=8.5)
ggsave("cumulative_lineage_counts_watchlist_mutations_b117", plot = plot_7, device = "pdf", width=11, height=8.5)


plot_list <- mget(ls(pattern = "^plot*"))
pdf("ncov_surveillance_all.pdf", width = 11, height = 8.5)
invisible(lapply(plot_list, print))
dev.off()


