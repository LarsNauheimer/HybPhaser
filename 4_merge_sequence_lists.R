#####################################################################
### Merging of sequence lists (non-phased with phased accessions) ###
#####################################################################

# load config
if (!(exists("config_file"))) {
	config_file <- "./config.txt"
}
source(config_file)

folder_seqlist_normal <- file.path(path_to_sequence_lists_normal)
folder_seqlist_phased <- file.path(path_to_sequence_lists_phased)
folder_output <- file.path(path_of_sequence_lists_output)

dir.create(folder_output, recursive = TRUE, showWarnings = FALSE)


seqlists_normal <- list.files(folder_seqlist_normal, full.names = TRUE)
seqlists_phased <- list.files(folder_seqlist_phased, full.names = TRUE)

# loci to include
if (file_with_loci_included == "") {
	loci_lists_normal <- gsub(".*/(.*)_c.*", "\\1", seqlists_normal)
	loci_lists_phased <- gsub(".*/(.*)_c.*", "\\1", seqlists_phased)
	loci2include <- unique(c(loci_lists_normal, loci_lists_phased))
} else {
	loci2include <- readLines(file_with_loci_included)
	loci2include <- gsub(" *$", "", loci2include)		# removing trailing spaces
}

# loci to exclude
if (file_with_loci_excluded != "") {
  loci2exclude <- readLines(file_with_loci_excluded)
  loci2exclude <- gsub(" *$", "", loci2exclude)			# removing trailing spaces
  loci2include <- loci2include[- which(loci2include %in% loci2exclude)]
}


# copy files from normal seq lists if in the list of loci to include
for (seqln in seqlists_normal) {
	if (gsub(".*/(.*)_c.*", "\\1", seqln) %in% loci2include) {
		file.copy(seqln, to = folder_output, overwrite = TRUE)
	}
}

seqlists_out_full <- list.files(folder_output, full.names = TRUE)
seqlists_out <- list.files(folder_output, full.names = FALSE)


for (seqlp in seqlists_phased) {
	loci_p <- gsub(".*/(.*)_c.*", "\\1", seqlp)
	if (loci_p %in% gsub(".*/(.*)_c.*", "\\1", seqlists_normal)) {
		file.append(seqlists_out_full[grep(paste0("\\b", loci_p, "\\b"),
			gsub("(.*)(_consensus.fasta|_contig.fasta)", "\\1", seqlists_out))], seqlp)
	} else if (include_phased_seqlists_when_non_phased_locus_absent != "no") {
		file.copy(seqlp, to = folder_output)
	}
}



# Samples - remove samples to be excluded from the sequence lists


# getting list included samples (either default from nameslists or from file_with_samples_included)
if (file_with_samples_included == "") {
	samples_include_normal <- readLines(path_to_namelist_normal)
	samples_include_phased <- readLines(path_to_namelist_phased)
	samples2include <- c(samples_include_normal, samples_include_phased)
} else {
	samples2include <- readLines(file_with_samples_included)
}
samples2include <- gsub(" *$", "", samples2include)		# removing trailing spaces


# exclude samples from file_with_samples_excluded
if (file_with_samples_excluded == "") {
	samples2exclude <- vector()
} else {
	samples2exclude <- readLines(file_with_samples_excluded)
	samples2exclude <- gsub(" *$", "", samples2exclude)	# removing trailing spaces
}


# removing unphased versions of phased samples, if that was selected
if (exchange_phased_with_not_phased_samples != "no") {
	phased_samples_base <- sub("_to_.*", "", samples_include_phased)
	notphased2remove <- samples2include[which(samples2include %in% phased_samples_base)]
} else {
	notphased2remove <- vector()
}

# remove excluded samples from the list of included samples
if (length(c(notphased2remove, samples2exclude)) > 0) {
	samples2include <- samples2include[- which(samples2include %in% c(notphased2remove, samples2exclude))]
}


# remove excluded samples from the sequences lists
for (seqlist_out in seqlists_out_full) {
	seqlist_out_raw <- readLines(seqlist_out)
	lines_with_samples_in <- which(seqlist_out_raw %in% paste0(">", samples2include))
	lines_to_keep <- sort(c(lines_with_samples_in, lines_with_samples_in + 1))
	seqlist_out_new <- seqlist_out_raw[lines_to_keep]
	conn <- file(seqlist_out)
	writeLines(seqlist_out_new, conn)
	close(conn)
}
