#########################################
### Preparation of the phasing script ###
#########################################

# load config
if (!(exists("config_file"))) {
	config_file <- "./config.txt"
}
source(config_file)


prep_phasing <- read.csv(csv_file_with_phasing_prep_info, header = TRUE)
prep_phasing[is.na(prep_phasing)] <- ""

path_to_reference_sequences <- file.path(reference_sequence_folder)
refseqs <- list.files(path_to_reference_sequences)
refseqs_fullpath <- list.files(path_to_reference_sequences, full.names = TRUE)
refseqs_samplenames <- gsub("(_intronerated)*_consensus.fasta|(_intronerated)*_contig.fasta", "", refseqs)
reads <- list.files(path_to_read_files_phasing, full.names = TRUE)

if (folder_for_phased_reads == "") {
	folder_for_phased_reads <- file.path(path_to_phasing_folder, "phased_reads")
}
dir.create(folder_for_phased_reads, recursive = TRUE, showWarnings = FALSE)

if (folder_for_phasing_stats == "") {
	folder_for_phasing_stats <- file.path(path_to_phasing_folder, "phasing_stats")
}
dir.create(folder_for_phasing_stats, recursive = TRUE, showWarnings = FALSE)


# prepare reference part of the command
ref_command <- vector()
for (i in seq_len(length(rownames(prep_phasing)))) {
	nrefs <- (length(which(prep_phasing[i, ] != "")) - 1) / 2
	ref_co <- vector()
	for (r in seq_len(nrefs)) {
		ref_file <- refseqs_fullpath[match(prep_phasing[i, r * 2], refseqs_samplenames)]
		ref_abb <- prep_phasing[i, (r * 2) + 1]
		ref_co[r] <- (paste0("ref_", ref_abb, "=", ref_file))
	}
	ref_command[i] <- paste(ref_co, collapse = " ")
}


# setting no of threads
if (no_of_threads_phasing == 0 || no_of_threads_phasing == "auto") {
	threadtext <- ""
} else {
	threadtext <- paste0(" threads=", no_of_threads_phasing)
}


# setting Java memory usage
if (java_memory_usage_phasing != "") {
	pXmx <- paste0(" -Xmx", java_memory_usage_phasing)
} else {
	pXmx <- ""
}

if (path_to_bbmap_executables == "") {
	bbsplit_sh <- "bbsplit.sh"
} else {
	bbsplit_sh <- file.path(path_to_bbmap_executables, "bbsplit.sh")
}


# preparing the whole command line
phasing_command <- vector()
if (read_type_4phasing == "paired-end") {
	for (i in seq_len(length(rownames(prep_phasing)))) {
		read_file_1 <- reads[match(paste0(prep_phasing[i, 1], ID_read_pair1), gsub(".*/", "", reads))]
		read_file_2 <- reads[match(paste0(prep_phasing[i, 1], ID_read_pair2), gsub(".*/", "", reads))]
		phasing_command[i] <- paste0(bbsplit_sh, " ambiguous=all ambiguous2=all", threadtext, " ",
			ref_command[i], " in=", read_file_1, " in2=", read_file_2, " basename=", folder_for_phased_reads,
			"/", prep_phasing[i, 1], "_to_%.fastq", " refstats=", folder_for_phasing_stats, "/",
			prep_phasing[i, 1], "_phasing-stats.txt", pXmx)
	}
} else if (read_type_4phasing == "single-end") {
	for (i in seq_len(length(rownames(prep_phasing)))) {
		read_file <- reads[grep(paste(prep_phasing[i, 1]), reads)]
		if (length(read_file) > 1) {
			stop(print("ERROR: multiple single end reads are selected!"), call. = FALSE)
		}
		phasing_command[i] <- paste0(bbsplit_sh, " ambiguous=all ambiguous2=all", threadtext, " ",
			ref_command[i], " in=", read_file, " basename=", folder_for_phased_reads, "/",
			prep_phasing[i, 1], "_to_%.fastq", " refstats=", folder_for_phasing_stats, "/",
			prep_phasing[i, 1], "_phasing-stats.txt", pXmx)
	}
} else {
	stop(print("ERROR: READ TYPE NOT SELECTED"), call. = FALSE)
}


#### generating phasing bash script
phasing_script_file <- file.path(path_to_phasing_folder, "run_bbsplit4phasing.sh")
write("#!/bin/bash", file = phasing_script_file, append = FALSE)
write(phasing_command, file = phasing_script_file, append = TRUE)
system(command = paste("chmod +x", phasing_script_file))
           
# run in R if selected
if (run_bash_script_in_R == "yes") {
	system(command = phasing_script_file)
}
