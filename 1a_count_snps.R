##################
### SNPs count ###
##################

# load config
if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)
library(seqinr)
library(stringr)

# generate folder

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
dir.create(output_Robjects, showWarnings = F, recursive = T)

#####################################################################################
### Counting polymorphic sites (masked as ambiguity codes in conseneus sequences) ###
#####################################################################################

# getting targets and sample names
if(targets_file_format == "AA"){
  targets <- read.fasta(fasta_file_with_targets, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
  } else if(targets_file_format == "DNA"){
    targets <- read.fasta(fasta_file_with_targets, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
  } else {
    print("Warning! Target file type not set properly. Should be 'DNA' or 'AA'!")
}
targets_name <- unique(gsub(".*-","",labels(targets)))

samples <- readLines(path_to_namelist)



if(intronerated_contig=="yes"){
  intronerated_name <- "intronerated" 
  intronerated_underscore <- "_"
} else {
  intronerated_name <- ""
  intronerated_underscore <- ""
  }

# function for counting ambiguities
seq_stats <- function(file){
  fasta <- read.fasta(file, as.string=TRUE, set.attributes = FALSE )
  seq <- gsub("N|[?]|-","",fasta[[1]])
  c(round(str_length(seq),0),(str_count(seq,"Y|K|R|S|M|y|k|r|s|m|w") + str_count(seq,"W|D|H|B|V|w|d|h|b|v")*2) )
}

# generate empty tables for SNPS and sequence length
tab_snps <- data.frame(loci=targets_name)
tab_length <- data.frame(loci=targets_name)

# fill tables with information on snps and sequence length for each sample and locus 
start_time <- Sys.time()
for(sample in samples){
  #sample <- samples[1]
  print(paste(sample," ", round(difftime(Sys.time(),start_time, units="secs"),0),"s", sep="") )
  
  tab_sample <- data.frame(targets=targets_name, seq_length=NA, ambis=NA, ambi_prop=NA)
  
  
  consensus_files <- list.files(file.path(path_to_output_folder,"01_data",sample, paste(intronerated_name,intronerated_underscore,"consensus", sep="")),full.names = T)
  
  for(consensus_file in consensus_files) {
    
    gene <- gsub("(_intronerated|).fasta","",gsub(".*/","",consensus_file))
    
    file.path(path_to_output_folder,"01_data",sample, "consensus", paste(gene,intronerated_underscore,intronerated_name,".fasta",sep=""))
    
    if(file.info(consensus_file)$size >1){
      stats <- seq_stats(consensus_file)
    } else {
      stats <- c(NA,NA)
    }  
    tab_sample$seq_length[match(gene,tab_sample$targets)] <- stats[1]
    tab_sample$ambis[match(gene,tab_sample$targets)] <- stats[2]  
  
  }
  
  tab_sample$ambi_prop <- tab_sample$ambis/tab_sample$seq_length  
  tab_snps[match(sample,samples)] <- tab_sample$ambi_prop
  tab_length[match(sample,samples)] <- tab_sample$seq_length

}


colnames(tab_snps) <- samples
rownames(tab_snps) <- targets_name
colnames(tab_length) <- samples
rownames(tab_length) <- targets_name

### Generate output tables and save data in Robjects


saveRDS(tab_snps,file=file.path(output_Robjects,"Table_SNPs.Rds"))
saveRDS(tab_length,file=file.path(output_Robjects,"Table_consensus_length.Rds"))
