###################################
### Collection of phasing stats ###
###################################

# load config
if (!(exists("config_file"))) {
	config_file <- "./config.txt"
}
source(config_file)


phasing_prep <- read.csv(csv_file_with_phasing_prep_info, header = TRUE)

if (folder_for_phasing_stats == "") {
	folder_for_phasing_stats <- file.path(path_to_phasing_folder, "phasing_stats/")
}
stats_files <- list.files(folder_for_phasing_stats, pattern = "*phasing-stats.txt")


# generate empty table (ref_names x samples)
used_reference_abbr <- unique(as.vector(t(phasing_prep[, 1 + 1:((length(phasing_prep) - 1) / 2) * 2])))
if (length(which(is.na(used_reference_abbr) | used_reference_abbr == "")) > 0) {
	used_reference_abbr <- used_reference_abbr[- which(is.na(used_reference_abbr) | used_reference_abbr == "")]
}


tab_phasing_stats <- matrix(nrow = length(stats_files), ncol = length(used_reference_abbr))
rownames(tab_phasing_stats) <- gsub("_phasing-stats.txt", "", stats_files)
colnames(tab_phasing_stats) <- used_reference_abbr


# fill table with results from the stats files
for (i in seq_len(length(stats_files))) {
	stats_file <- stats_files[i]
	p_unamb <- read.delim(file.path(folder_for_phasing_stats, stats_file))[, 1:2]
  for (j in seq_len(length(p_unamb[, 1]))) {
		tab_phasing_stats[i, match(p_unamb[j, 1], colnames(tab_phasing_stats))] <- round(p_unamb[j, 2] / 100, 5)
  }
}


# save output
write.csv(tab_phasing_stats,
	file = file.path(path_to_phasing_folder, "Table_phasing_stats.csv"), na = "")
