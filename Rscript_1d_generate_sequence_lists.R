
#############################################
#############################################

if(contig=="supercontig"){
  contig_name <- "_supercontig"
  contig_consensus_name <- "supercontig-"
  output_base <- file.path(path_to_hybpiper_results,paste("_HybPhaser_supercontig",sep=""))
} else {
  contig_name <- ""
  contig_consensus_name <- ""
  output_base <- file.path(path_to_hybpiper_results,paste("_HybPhaser",sep=""))
}

targets <- read.fasta(fasta_file_with_targets, as.string=TRUE, set.attributes = FALSE)
targets_name <- unique(gsub(".*-","",labels(targets)))
samples <- readLines(txt_file_with_list_of_accessions)

#outloci_para_all <- readLines(file.path(output_base,"cleaning1_loci_with_high_divergence_across_all_samples.txt"))
outsamples_missing <- readRDS(file=file.path(output_Robjects,"outsamples_missing.Rds"))
outloci_missing <- readRDS(file=file.path(output_Robjects,"outloci_missing.Rds"))
outloci_para_all <- readRDS(file=file.path(output_Robjects,"outloci_para_all.Rds"))
outloci_para_each <- readRDS(file=file.path(output_Robjects,"outloci_para_each.Rds"))
tab_snps_cl2b <- readRDS(file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))

output_sequences <- file.path(path_for_HybPhaser_output,"sequence_lists/")


#############################################


folder4seq_consensus_loci_raw <- file.path(output_sequences,"consensus_loci_raw")
folder4seq_consensus_loci_clean <- file.path(output_sequences,"consensus_loci_clean")
folder4seq_consensus_samples_raw <- file.path(output_sequences,"consensus_samples_raw")
folder4seq_consensus_samples_clean <- file.path(output_sequences,"consensus_samples_clean")
folder4seq_contig_loci_raw <- file.path(output_sequences,"contig_loci_raw")
folder4seq_contig_loci_clean <- file.path(output_sequences,"contig_loci_clean")
folder4seq_contig_samples_raw <- file.path(output_sequences,"contig_samples_raw")
folder4seq_contig_samples_clean <- file.path(output_sequences,"contig_samples_clean")

unlink(output_sequences,recursive = T) # delete directory, if it existed in order to prevent errors
dir.create(output_sequences, showWarnings = T)
dir.create(folder4seq_consensus_loci_raw, showWarnings = T)
dir.create(folder4seq_consensus_loci_clean, showWarnings = T)
dir.create(folder4seq_consensus_samples_raw, showWarnings = T)
dir.create(folder4seq_consensus_samples_clean, showWarnings = T)
dir.create(folder4seq_contig_loci_raw, showWarnings = T)
dir.create(folder4seq_contig_loci_clean, showWarnings = T)
dir.create(folder4seq_contig_samples_raw, showWarnings = T)
dir.create(folder4seq_contig_samples_clean, showWarnings = T)




#########################################################################
### concatenate consensus files across all samples to lists per locus ###
#########################################################################
# this extracts sequences from all subfolders in the HybPiper folder and collates them into one file per locus

# HybPhaser
if(contig!="supercontig"){
  for(locus in targets_name){
    command_cat_consensus <- paste("cat",file.path(path_to_hybpiper_results,"*",locus,"*/sequences/remapping/",paste(locus,"_consensus.fasta",sep="")),">",file.path(folder4seq_consensus_loci_raw,paste(locus,"_consensus_raw.fasta",sep="")))
    system(command_cat_consensus)
    command_cat_contig <- paste("cat",file.path(path_to_hybpiper_results,"*",locus,"*/sequences/FNA/",paste(locus,".FNA",sep="")),">",file.path(folder4seq_contig_loci_raw,paste(locus,"_contig_raw.fasta",sep="")))
    system(command_cat_contig)
    command_remove_locus_in_seqnames <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci_raw,paste(locus,"_consensus_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames)
    command_remove_locus_in_seqnames <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci_raw,paste(locus,"_contig_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames)
  } 
  } else {
    for(locus in targets_name){
      command_cat_consensus <- paste("cat",file.path(path_to_hybpiper_results,"*",locus,"*/sequences/remapping/",paste(locus,"_supercontig-consensus.fasta",sep="")),">",file.path(folder4seq_consensus_loci_raw,paste(locus,"_supercontig-consensus_raw.fasta",sep="")))
      system(command_cat_consensus)
      command_cat_contig <- paste("cat",file.path(path_to_hybpiper_results,"*",locus,"*/sequences/intron/",paste(locus,"_supercontig.fasta",sep="")),">",file.path(folder4seq_contig_loci_raw,paste(locus,"_supercontig-contig_raw.fasta",sep="")))
      system(command_cat_contig)
      command_remove_locus_in_seqnames <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci_raw,paste(locus,"_supercontig-consensus_raw.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames) 
      command_remove_locus_in_seqnames <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci_raw,paste(locus,"_supercontig-contig_raw.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames)
    }
}  
  
for(file in list.files(folder4seq_contig_loci_raw, full.names = T)){
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}



########################################################################
### remove loci from dtaset optimization (missing data and paralogs) ###
########################################################################

## copy all sequence lists to new folder before removing relevant sequences

# consensus
file.remove(list.files(folder4seq_consensus_loci_clean,full.names = T)) # clear folder
file.copy(from=list.files(folder4seq_consensus_loci_raw,full.names = T), to = folder4seq_consensus_loci_clean) # copy raw files
system(paste("rename 's/_raw.fasta/.fasta/' ",folder4seq_consensus_loci_clean,"/*raw.fasta -f", sep=""))  # rename files (removing "_raw")

# contig
file.remove(list.files(folder4seq_contig_loci_clean,full.names = T)) # clear folder
file.copy(from=list.files(folder4seq_contig_loci_raw,full.names = T), to = folder4seq_contig_loci_clean) # copy raw files
system(paste("rename 's/_raw.fasta/.fasta/' ",folder4seq_contig_loci_clean,"/*raw.fasta -f", sep=""))  # rename files (removing "_raw")


## remove all outlier loci

files_consensus_clean <- list.files(path = folder4seq_consensus_loci_clean, full.names = T )
files_contig_clean <- list.files(path = folder4seq_contig_loci_clean, full.names = T )

loci_to_remove <- c(failed_loci, outloci_missing, outloci_para_all)

if(length(loci_to_remove)==0){
  files_to_remove_consensus=""
  files_to_remove_contig=""
} else {
  files_to_remove_consensus <- files_consensus_clean[grep(paste(c(names(loci_to_remove)),collapse="|"),files_consensus_clean)]
  files_to_remove_contig <- files_contig_clean[grep(paste(c(names(loci_to_remove)),collapse="|"),files_contig_clean)]
}

file.remove(files_to_remove_consensus)
file.remove(files_to_remove_contig)



## remove sequences in locus files from paralogs for each sample


for(locus in rownames(tab_snps_cl2b)){

  samples_to_remove <- vector()
  
  if(length(failed_samples) > 0 ){
    samples_to_remove <- names(failed_samples)
  } 
  
  if(length(outsamples_missing) > 0 ){
    samples_to_remove <-  c(samples_to_remove, names(outsamples_missing))
  } 
  
  if(length(outsamples_recovered_seq_length) > 0 ){
    samples_to_remove <- c(samples_to_remove, names(outsamples_recovered_seq_length))
  } 
  
  if(length(grep(locus,outloci_para_each)) >0 ){
    samples_to_remove <- c(samples_to_remove, names(outloci_para_each[grep(locus,outloci_para_each)]))
  }
  
  
  if(length(samples_to_remove) !=0 ){
    
    # consensus
    locus_consensus_raw <- readLines(file.path(folder4seq_consensus_loci_clean,paste(locus,"_",contig_consensus_name,"consensus.fasta",sep="")))
    lines_with_samplename <- grep(paste(samples_to_remove,collapse="|"),locus_consensus_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus_raw[-lines_to_remove]
      conn <- file(file.path(folder4seq_consensus_loci_clean,paste(locus,"_",contig_consensus_name,"consensus.fasta",sep="")))
      writeLines(locus_file_red, conn)
      close(conn)
    }
    
    # contig
    locus_contig_raw <- readLines(file.path(folder4seq_contig_loci_clean,paste(locus,contig_name,"-contig.fasta",sep="")))
    lines_with_samplename <- grep(paste(samples_to_remove,collapse="|"),locus_contig_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig_raw[-lines_to_remove]
      conn <- file(file.path(folder4seq_contig_loci_clean,paste(locus,contig_name,"-contig.fasta",sep="")))
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
    
  }
} 



#########################################################################
### concatenate consensus files across all loci to lists per sample   ###
#########################################################################

samples_in <- samples
if(length(c(failed_samples, outsamples_missing, outsamples_recovered_seq_length)) != 0){
  samples_in <- samples[-which(samples %in% unique(c(names(failed_samples), names(outsamples_missing) , names(outsamples_recovered_seq_length))))]
}


if(contig=="supercontig"){
  for(sample in samples){
    command_cat_loci_consensus <- paste("cat",file.path(path_to_hybpiper_results,sample,"/*",sample,"/sequences/remapping/*_supercontig-consensus.fasta"),">",file.path(folder4seq_consensus_samples_raw,paste(sample,"_supercontig-consensus_raw.fasta",sep="")))
    system(command_cat_loci_consensus)
    command_cat_loci_contig <- paste("cat",file.path(path_to_hybpiper_results,sample,"/*",sample,"/sequences/intron/*_supercontig.fasta"),">",file.path(folder4seq_contig_samples_raw,paste(sample,"_contig_supercontig_raw.fasta",sep="")))
    system(command_cat_loci_contig)
  } 
} else {
  for(sample in samples){
    # for consensus files
    command_cat_loci_consensus <- paste("cat",file.path(path_to_hybpiper_results,sample,"/*",sample,"/sequences/remapping/*_consensus.fasta"),">",file.path(folder4seq_consensus_samples_raw,paste(sample,"_consensus_raw.fasta",sep="")))
    system(command_cat_loci_consensus)
    #for contig files (they cont contain the gene name , which makes it more complicated)
    dir.create(file.path(folder4seq_contig_samples_raw,"tmp"))
    system(paste("find ", file.path(path_to_hybpiper_results,sample)," -type f -wholename '*/FNA/*.FNA' -exec cp {} ", file.path(folder4seq_contig_samples_raw,"tmp")," \\;"))
    tmp_files <- list.files(file.path(folder4seq_contig_samples_raw,"tmp"), full.names = T)
    for(filename in tmp_files){
      gene_name <- gsub(".*/(.*).FNA","\\1",filename)
      system(paste("sed -i 's/",sample, "/",gene_name,"/g' ", filename, sep=""))
    }
    system(paste("cat",file.path(folder4seq_contig_samples_raw,"tmp/*.*"),">",file.path(folder4seq_contig_samples_raw,paste(sample,"_contig_raw.fasta",sep=""))))
    unlink(file.path(folder4seq_contig_samples_raw,"tmp"), recursive = T)
    
  }
}  





dir.create(file.path(folder4seq_contig_samples_raw,"tmp"))
system(cat("find ", file.path(path_to_hybpiper_results,sample)," -type f -wholename '*/FNA/*.FNA' -exec cp {} ", file.path(folder4seq_contig_samples_raw,"tmp")," \\;"))
tmp_files <- list.files(file.path(folder4seq_contig_samples_raw,"tmp"), full.names = T)
for(filename in tmp_files){
  gene_name <- gsub(".*/(.*).FNA","\\1",filename)
  system(paste("sed -i 's/",sample, "/",gene_name,"/g' ", filename, sep=""))
}
system(paste("cat",file.path(folder4seq_contig_samples_raw,"tmp/*.*"),">",file.path(folder4seq_contig_samples_raw,paste(sample,"_contig_raw.fasta",sep=""))))
unlink(file.path(folder4seq_contig_samples_raw,"tmp"), recursive = T)





# changing interelaved HybPiper files to non-interleaves fasta files
for(file in list.files(folder4seq_contig_samples_raw, full.names = T)){
    lines <- readLines(file)
    conx <- file(file)
    writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
    close(conx)
}




###########################################################################
## removing loci from sample lists
###########################################################################

## copy all sequence lists to new folder before removing parts of it

# consensus
file.remove(list.files(folder4seq_consensus_samples_clean,full.names = T)) # clear folder
file.copy(from=list.files(folder4seq_consensus_samples_raw,full.names = T), to = folder4seq_consensus_samples_clean)
system(paste("rename 's/_raw.fasta/.fasta/' ",folder4seq_consensus_samples_clean,"/*raw.fasta -f", sep=""))

#contig
file.remove(list.files(folder4seq_contig_samples_clean,full.names = T)) # clear folder
file.copy(from=list.files(folder4seq_contig_samples_raw,full.names = T), to = folder4seq_contig_samples_clean)
system(paste("rename 's/_raw.fasta/.fasta/' ",folder4seq_contig_samples_clean,"/*raw.fasta -f", sep=""))


## remove all outlier loci (empty and high allele divergence) (cleaning step 1 and 2)
sample_files_consensus_clean <- list.files(path = folder4seq_consensus_samples_clean, full.names = T )
sample_files_contig_clean <- list.files(path = folder4seq_contig_samples_clean, full.names = T )

for(sample in samples_in){
  
  loci_to_remove <- vector()
  
  if(length(failed_loci) > 0){
    loci_to_remove <- names(failed_loci)
  }
  
  if(length(outloci_missing) > 0){
    loci_to_remove <- c(loci_to_remove, names(outloci_missing))
  }
  
  if(length(outloci_para_all) > 0){
    loci_to_remove <- c(loci_to_remove, names(outloci_para_all))
  }
  
  if(length(grep(sample,names(outloci_para_each))) > 0 ){
    loci_to_remove <- c(loci_to_remove, outloci_para_each[[which(names(outloci_para_each) %in% sample)]])
  }
  

  if(length(samples_to_remove)!=0){
    
    
    #consensus
    consensus_file2clean <- file.path(folder4seq_consensus_samples_clean, paste(sample,"_",contig_consensus_name,"consensus.fasta",sep=""))
    samples_consensus_raw <- readLines(consensus_file2clean)
    lines_with_lociname <- grep(paste(loci_to_remove,collapse="|"),samples_consensus_raw)
    if(length(lines_with_lociname) !=0){
      lines_to_remove <- c(lines_with_lociname,lines_with_lociname+1)
      sample_file_consensus_red <- samples_consensus_raw[-lines_to_remove]
      conn <- file(consensus_file2clean)
      writeLines(sample_file_consensus_red, conn)
      close(conn)
    }
    
    #contig
    contig_file2clean <- file.path(folder4seq_contig_samples_clean, paste(sample,contig_name,"_contig.fasta",sep=""))
    samples_contig_raw <- readLines(contig_file2clean)
    lines_with_lociname <- grep(paste(loci_to_remove,collapse="|"),samples_contig_raw)
    if(length(lines_with_lociname) !=0){
      lines_to_remove <- c(lines_with_lociname,lines_with_lociname+1)
      sample_file_contig_red <- samples_contig_raw[-lines_to_remove]
      conn <- file(contig_file2clean)
      writeLines(sample_file_contig_red, conn)
      close(conn)
    }
  }
} 

