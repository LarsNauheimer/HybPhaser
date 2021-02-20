################################
# Assigning folder structure ###
################################


dir.create(path_for_HybPhaser_output, showWarnings = F)

output_cleaning <- file.path(path_for_HybPhaser_output,"dataset_optimization/")
output_assess <- file.path(path_for_HybPhaser_output,"heterozygosity_assessment/")
output_Robjects <- file.path(path_for_HybPhaser_output,"R_objects/")
output_sequences <- file.path(path_for_HybPhaser_output,"sequence_lists/")


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

samples <- readLines(txt_file_with_list_of_accessions)


if(contig=="supercontig"){consensus_name <- "_supercontig-consensus" 
} else {consensus_name <- "_consensus"}

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
  
  good_genes <- vector() 
  if(file.exists(file.path(path_to_hybpiper_results,sample,"genes_with_seqs.txt"))){good_genes <- sort(sub("\t.*","",readLines(file.path(path_to_hybpiper_results,sample,"genes_with_seqs.txt"))))}
  
  for(gene in good_genes) {
    
    consensus_file <- file.path(path_to_hybpiper_results,sample,gene,sample,"sequences/remapping",paste(gene,consensus_name,".fasta",sep=""))
    
    if(file.exists(consensus_file)){
      if(file.info(consensus_file)$size !=0){
        stats <- seq_stats(consensus_file)
      } else {
        stats <- c(NA,NA)
      }  
    } else {
      stats <- c(NA,NA)
    } 
    
    tab_sample$seq_length[grep(gene,tab_sample$targets)] <- stats[1]
    tab_sample$ambis[grep(gene,tab_sample$targets)] <- stats[2]  
    

  }
  
  tab_sample$ambi_prop <- tab_sample$ambis/tab_sample$seq_length  
  tab_snps[grep(sample,samples)+1] <- tab_sample$ambi_prop
  tab_length[grep(sample,samples)+1] <- tab_sample$seq_length

}

colnames(tab_snps) <- c("loci",samples)
colnames(tab_length) <- c("loci",samples)


### Generate output tables and save data in Robjects
dir.create(output_assess, showWarnings = F)

write.csv(tab_snps, file = file.path(output_assess,"Table_SNPs_raw.csv"), row.names = tab_snps$loci)
write.csv(tab_length, file = file.path(output_assess,"Table_consensus_length.csv"), row.names = tab_snps$loci)

dir.create(output_Robjects, showWarnings = F)

saveRDS(tab_snps,file=file.path(output_Robjects,"Table_SNPs_raw.Rds"))
saveRDS(tab_length,file=file.path(output_Robjects,"Table_consensus_length.Rds"))


