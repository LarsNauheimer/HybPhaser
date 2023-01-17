###########################################################
### Preparation of BBSplit script for Clade association ###
###########################################################

# load config
if (!(exists("config_file"))) {
	config_file <- "./config.txt"
}
source(config_file)


if (read_type_cladeassociation == "paired-end") {
	read_files <- list.files(path_to_read_files_cladeassociation, full.names = FALSE,
		pattern = paste0("*", ID_read_pair1, ".*"), include.dirs = FALSE)
} else if (read_type_cladeassociation == "single-end") {
	read_files <- list.files(path_to_read_files_cladeassociation, full.names = FALSE,
		pattern = "*.*", include.dirs = FALSE)
  	ID_read_pair1 <- ""
} else {
	stop(print("ERROR: READ TYPE NOT SELECTED"), call. = FALSE)
}

if (file_with_samples_included != "" && file_with_samples_included != "none") {
	samples_in <- readLines(file_with_samples_included)
	samples_in <- gsub("\\s*$", "", samples_in)
	read_files_noID <- gsub(ID_read_pair1, "", read_files)
	read_files_noend <- gsub("[.]fast.*", "", read_files_noID)
	read_files <- read_files[match(samples_in, read_files_noend)]
}


read_files <- na.exclude(read_files)

folder_bbsplit_stats <- file.path(path_to_clade_association_folder, "bbsplit_stats/")
dir.create(folder_bbsplit_stats, showWarnings = FALSE)

sample_sequences <- list.files(path_to_reference_sequences)

ref_samples <- read.csv(csv_file_with_clade_reference_names, header = TRUE)
colnames(ref_samples)[1] <- "samples"
colnames(ref_samples)[2] <- "abb"

ref_samples_files <- sample_sequences[match(ref_samples$samples,
	gsub("(_intronerated|_consensus|_contig).*", "", sample_sequences))]

if (length(ref_samples$abb) > 0) {
	ref_command <- paste(paste0("ref_", ref_samples[, 2], "=", path_to_reference_sequences, "/", ref_samples_files),
		collapse = " ")
} else {
	ref_command <- paste("ref=", paste0(path_to_reference_sequences, ref_samples_files, collapse = ","))
}


# setting no of threads
if (no_of_threads_clade_association == 0 || no_of_threads_clade_association == "auto") {
	threadtext <- ""
} else {
	threadtext <- paste0(" threads=", no_of_threads_clade_association)
}


# setting Java memory usage
if (java_memory_usage_clade_association != "") {
	caXmx <- paste0(" -Xmx", java_memory_usage_clade_association)
} else {
	caXmx <- ""
}

if (path_to_bbmap == "") {
	bbsplit_sh <- "bbsplit.sh"
} else {
	bbsplit_sh <- file.path(path_to_bbmap, "bbsplit.sh")
}



write("#!/bin/bash", file = file.path(path_to_clade_association_folder, "run_bbsplit4clade_association.sh"),
	append = FALSE)

for (i in seq_len(length(read_files))) {
	if (read_type_cladeassociation == "paired-end") {
		bbsplit_command <- paste0(bbsplit_sh, " ", ref_command, " in=",
			file.path(path_to_read_files_cladeassociation, read_files[i]), " in2=",
			sub(ID_read_pair1, ID_read_pair2, file.path(path_to_read_files_cladeassociation, read_files[i])),
			threadtext, " ambiguous=random ambiguous2=all refstats=", folder_bbsplit_stats,
			gsub(ID_read_pair1, "", read_files[i]), "_bbsplit-stats.txt", caXmx)
	} else {
		bbsplit_command <- paste0(bbsplit_sh, " ", ref_command, " in=",
			file.path(path_to_read_files_cladeassociation, read_files[i]), threadtext,
			" ambiguous=random ambiguous2=all refstats=", folder_bbsplit_stats, read_files[i],
			"_bbsplit-stats.txt", caXmx)
	}
	write(bbsplit_command, file = file.path(path_to_clade_association_folder, "run_bbsplit4clade_association.sh"),
		append = TRUE)
}

system(command = paste("chmod +x", file.path(path_to_clade_association_folder, "run_bbsplit4clade_association.sh")))

if (run_clade_association_mapping_in_R == "yes") {
	system(command = file.path(path_to_clade_association_folder, "run_bbsplit4clade_association.sh"))
}
