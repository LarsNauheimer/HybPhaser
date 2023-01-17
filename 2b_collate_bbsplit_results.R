##############################################################
### Collection of info in log files from clade association ###
##############################################################

# load config
if (!(exists("config_file"))) {
	config_file <- "./config.txt"
}
source(config_file)

# get names of stats files and from it get the name of samples
folder_bbsplit_stats <- file.path(path_to_clade_association_folder, "bbsplit_stats")
stats_files <- list.files(folder_bbsplit_stats, "*_bbsplit-stats.txt")
# remove everything after the first . in the name
sample <- gsub("(_mapped.*|.fasta.*|_bbsplit-stats.*)", "", stats_files)

ref_samples <- read.csv(csv_file_with_clade_reference_names, header = TRUE)
colnames(ref_samples)[1] <- "sample"
colnames(ref_samples)[2] <- "abbreviation"

# set reference names (full names or abbreviations)
if ("abbreviation" %in% colnames(ref_samples)) {
	ref_names <- ref_samples$abbreviation
} else {
	ref_names <- ref_samples$sample
}


# generate empty table (ref_names x samples)
tab_clade_assoc <- matrix(nrow = length(ref_names), ncol = length(sample))
rownames(tab_clade_assoc) <- ref_names
colnames(tab_clade_assoc) <- sample


for (i in seq_len(length(stats_files))) {
	stats_file <- stats_files[i]
	p_unamb <- read.delim(file.path(folder_bbsplit_stats, stats_file))[, 1:2]
	tab_clade_assoc[, i] <- p_unamb$X.unambiguousReads[match(rownames(tab_clade_assoc), p_unamb$X.name)]
}


# The proportions of reads matching unambigously to a reference are sensitive to the genetic distance of references.
# E.g. two clade references that are closer related to each other have lower proportions of reads matching
# to the reference because more reads match ambiguously between those references.
# One way to get more comparable data is normalization of proportions, by dividing the proportion of reads matched
# to a reference by the results from the sample that was used as reference.
# This assumes that the reads of the reference sequence matching to the same reference is the maximum possible and
# therefore set to 1.
# However, this can deviate due to differences in sequencing quality, read coverage and other factors.


tab_clade_assoc_norm <- tab_clade_assoc

for (i in seq_len(length(rownames(tab_clade_assoc)))) {
	if (length(grep(paste0("\\b", ref_samples$sample[i], "\\b"), colnames(tab_clade_assoc))) > 0) {
		numerator <- tab_clade_assoc[i, ]
		denominator <- tab_clade_assoc[i, grep(paste0("\\b", ref_samples$sample[i], "\\b"), colnames(tab_clade_assoc))]
		tab_clade_assoc_norm[i, ] <- numerator / denominator
  }
}

match("3", c("1", "2", "22", "3"))

ttab_clade_assoc <- as.data.frame(t(tab_clade_assoc))
ttab_clade_assoc$sample <- rownames(ttab_clade_assoc)

ttab_clade_assoc_norm <- as.data.frame(t(tab_clade_assoc_norm))
ttab_clade_assoc_norm$sample <- rownames(ttab_clade_assoc_norm)


# for better assessment, the summary table is added to the clade assessment results
summary_het_ad_clean <- readRDS(file.path(path_to_output_folder, "00_R_objects",
	name_for_dataset_optimization_subset, "Summary_table.Rds"))

tab_hetad_cladeasso <- merge(summary_het_ad_clean, ttab_clade_assoc, by = "sample", incomparables = FALSE)

tab_hetad_cladeasso_norm <- merge(summary_het_ad_clean, ttab_clade_assoc_norm, by = "sample")

# save output
write.csv(t(tab_clade_assoc),
	file = file.path(path_to_clade_association_folder, "Table_clade_association.csv"))
write.csv(t(tab_clade_assoc_norm),
	file = file.path(path_to_clade_association_folder, "Table_clade_association_normalised.csv"))
write.csv(tab_hetad_cladeasso,
	file = file.path(path_to_clade_association_folder, "Table_clade_association_and_summary_table.csv"))
write.csv(tab_hetad_cladeasso_norm,
	file = file.path(path_to_clade_association_folder, "Table_clade_association_normalized_and_summary_table.csv"))
