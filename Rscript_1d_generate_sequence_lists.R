
#############################################
#############################################

if(contig=="supercontig"){
  contig_name <- "supercontig-"
} else {
  contig_name <- ""
}

targets <- read.fasta(fasta_file_with_targets, as.string=TRUE, set.attributes = FALSE)
targets_name <- unique(gsub(".*-","",labels(targets)))
samples <- readLines(txt_file_with_list_of_accessions)

outsamples_missing <- readRDS(file=file.path(output_Robjects,"outsamples_missing.Rds"))
outloci_missing <- readRDS(file=file.path(output_Robjects,"outloci_missing.Rds"))
outloci_para_all <- readRDS(file=file.path(output_Robjects,"outloci_para_all.Rds"))
outloci_para_each <- readRDS(file=file.path(output_Robjects,"outloci_para_each.Rds"))
tab_snps_cl2b <- readRDS(file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))

output_sequences <- file.path(path_for_HybPhaser_output,"sequence_lists")


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
    command_cat_contig    <- paste("cat",file.path(path_to_hybpiper_results,"*",locus,"*/sequences/FNA/",paste(locus,".FNA",sep="")),">",file.path(folder4seq_contig_loci_raw,paste(locus,"_contig_raw.fasta",sep="")))
    system(command_cat_consensus)
    system(command_cat_contig)
    command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci_raw,paste(locus,"_consensus_raw.fasta",sep="")), sep=""))
    command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci_raw,paste(locus,"_contig_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames_consensus)
    system(command_remove_locus_in_seqnames_contig)
  } 
} else {
  for(locus in targets_name){
    command_cat_consensus <- paste("cat",file.path(path_to_hybpiper_results,"*",locus,"*/sequences/remapping/",paste(locus,"_supercontig-consensus.fasta",sep="")),">",file.path(folder4seq_consensus_loci_raw,paste(locus,"_supercontig-consensus_raw.fasta",sep="")))
    system(command_cat_consensus)
    command_cat_contig <- paste("cat",file.path(path_to_hybpiper_results,"*",locus,"*/sequences/intron/",paste(locus,"_supercontig.fasta",sep="")),">",file.path(folder4seq_contig_loci_raw,paste(locus,"_supercontig-contig_raw.fasta",sep="")))
    system(command_cat_contig)
    command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci_raw,paste(locus,"_supercontig-consensus_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames_consensus) 
    command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci_raw,paste(locus,"_supercontig-contig_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames_contig)
  }
}  

for(file in list.files(folder4seq_contig_loci_raw, full.names = T)){
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}


# check whether accessions are in the hybpiper folder but not in the samples list.
# if the sample list is smaller some sequences have to be removed from the loci lists

hybpiper_result_dirs <- list.dirs(path_to_hybpiper_results, full.names = FALSE, recursive = FALSE)
dirs_not_in_sample_list <- hybpiper_result_dirs[which(!(hybpiper_result_dirs %in% samples))]

if(length(dirs_not_in_sample_list) !=0 ){
  # consensus
  for(raw_consensus_file in list.files(folder4seq_consensus_loci_raw, full.names = TRUE)){
    locus_consensus_raw <- readLines(raw_consensus_file)
    lines_with_samplename <- grep(paste(dirs_not_in_sample_list,collapse="|"),locus_consensus_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus_raw[-lines_to_remove]
      conn <- file(raw_consensus_file)
      writeLines(locus_file_red, conn)
      close(conn)
    }
  }
  # contig
  for(raw_contig_file in list.files(folder4seq_contig_loci_raw, full.names = TRUE)){
    locus_contig_raw <- readLines(raw_contig_file)
    lines_with_samplename <- grep(paste(dirs_not_in_sample_list,collapse="|"),locus_contig_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig_raw[-lines_to_remove]
      conn <- file(raw_contig_file)
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
  }
}



### remove loci from dataset optimization (missing data and paralogs)
######################################################################

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

loci_files_consensus_clean <- list.files(path = folder4seq_consensus_loci_clean, full.names = T )
loci_files_contig_clean <- list.files(path = folder4seq_contig_loci_clean, full.names = T )

loci_to_remove <- c(failed_loci, outloci_missing, names(outloci_para_all))

if(length(loci_to_remove)==0){
  loci_files_to_remove_consensus=""
  loci_files_to_remove_contig=""
} else {
  loci_files_to_remove_consensus <- loci_files_consensus_clean[grep(paste(c(loci_to_remove),collapse="|"),loci_files_consensus_clean)]
  loci_files_to_remove_contig <- loci_files_contig_clean[grep(paste(c(loci_to_remove),collapse="|"),loci_files_contig_clean)]
}

file.remove(loci_files_to_remove_consensus)
file.remove(loci_files_to_remove_contig)



## remove sequences in locus files from paralogs for each sample


for(locus in rownames(tab_snps_cl2b)){
  
  samples_to_remove <- vector()
  
  if(length(failed_samples) > 0 ){
    samples_to_remove <- names(failed_samples)
  } 
  
  if(length(outsamples_missing) > 0 ){
    samples_to_remove <-  unique(c(samples_to_remove, outsamples_missing))
  } 
  
  if(length(grep(locus,outloci_para_each)) >0 ){
    samples_to_remove <- c(samples_to_remove, names(outloci_para_each[grep(locus,outloci_para_each)]))
  }
  
  
  if(length(samples_to_remove) !=0 ){
    
    # consensus
    locus_consensus_raw <- readLines(file.path(folder4seq_consensus_loci_clean,paste(locus,"_",contig_name,"consensus.fasta",sep="")))
    lines_with_samplename <- grep(paste(samples_to_remove,collapse="|"),locus_consensus_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus_raw[-lines_to_remove]
      conn <- file(file.path(folder4seq_consensus_loci_clean,paste(locus,"_",contig_name,"consensus.fasta",sep="")))
      writeLines(locus_file_red, conn)
      close(conn)
    }
    
    # contig
    locus_contig_raw <- readLines(file.path(folder4seq_contig_loci_clean,paste(locus,"_",contig_name,"contig.fasta",sep="")))
    lines_with_samplename <- grep(paste(samples_to_remove,collapse="|"),locus_contig_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig_raw[-lines_to_remove]
      conn <- file(file.path(folder4seq_contig_loci_clean,paste(locus,"_",contig_name,"contig.fasta",sep="")))
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
    
  }
} 



#########################################################################
### concatenate consensus files across all loci to lists per sample   ###
#########################################################################

samples <- readLines(txt_file_with_list_of_accessions)

# remove failed samples from list
if(length(failed_samples) != 0){
  samples <- samples[-which(samples %in% names(failed_samples))]
}

# collect all sequences

if(contig=="supercontig"){
  for(sample in samples){
    command_cat_loci_consensus <- paste("cat",file.path(path_to_hybpiper_results,sample,"/*",sample,"/sequences/remapping/*_supercontig-consensus.fasta"),">",file.path(folder4seq_consensus_samples_raw,paste(sample,"_supercontig-consensus_raw.fasta",sep="")))
    system(command_cat_loci_consensus)
    command_cat_loci_contig <- paste("cat",file.path(path_to_hybpiper_results,sample,"/*",sample,"/sequences/intron/*_supercontig.fasta"),">",file.path(folder4seq_contig_samples_raw,paste(sample,"_supercontig-contig_raw.fasta",sep="")))
    system(command_cat_loci_contig)
  } 
} else {
  for(sample in samples){
    # for consensus files
    command_cat_loci_consensus <- paste("cat",file.path(path_to_hybpiper_results,sample,"/*",sample,"/sequences/remapping/*_consensus.fasta"),">",file.path(folder4seq_consensus_samples_raw,paste(sample,"_consensus_raw.fasta",sep="")))
    system(command_cat_loci_consensus)
    #for contig files (they cont contain the gene name , which makes it more complicated)
    if(write_gene_names_in_contig_sample_seqlist!="yes"){
      command_cat_loci_contig <- paste("cat",file.path(path_to_hybpiper_results,sample,"/*",sample,"/sequences/FNA/*.FNA"),">",file.path(folder4seq_contig_samples_raw,paste(sample,"_contig_raw.fasta",sep="")))
      system(command_cat_loci_contig)
    }else {
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
}  


# changing interelaved HybPiper files to non-interleaves fasta files
for(file in list.files(folder4seq_contig_samples_raw, full.names = T)){
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}





# remove samples (missing data, paralogs for all)
###################################################

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


if(length(outsamples_missing)==0){
  samples_files_to_remove_consensus=""
  samples_files_to_remove_contig=""
} else {
  samples_files_to_remove_consensus <- sample_files_consensus_clean[grep(paste(c(outsamples_missing),collapse="|"),sample_files_consensus_clean)]
  samples_files_to_remove_contig <- sample_files_contig_clean[grep(paste(c(outsamples_missing),collapse="|"),sample_files_contig_clean)]
}

file.remove(samples_files_to_remove_consensus)
file.remove(samples_files_to_remove_contig)


samples_in <- samples
if(length(outsamples_missing) != 0){
  samples_in <- samples_in[-which(samples %in% outsamples_missing)]
}

# remove outlier loci per sample in sample lists
for(sample in samples_in){
  
  loci_to_remove <- vector()
  
  if(length(failed_loci) > 0){
    loci_to_remove <- failed_loci
  }
  
  if(length(outloci_missing) > 0){
    loci_to_remove <- c(loci_to_remove, outloci_missing)
  }
  
  if(length(outloci_para_all) > 0){
    loci_to_remove <- c(loci_to_remove, names(outloci_para_all))
  }
  
  if(length(grep(sample,names(outloci_para_each))) > 0 ){
    loci_to_remove <- c(loci_to_remove, names(outloci_para_each[[which(names(outloci_para_each) %in% sample)]]))
  }
  
  
  if(length(loci_to_remove)!=0){
    
    
    #consensus
    consensus_file2clean <- file.path(folder4seq_consensus_samples_clean, paste(sample,"_",contig_name,"consensus.fasta",sep=""))
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
    contig_file2clean <- file.path(folder4seq_contig_samples_clean, paste(sample,"_",contig_name,"contig.fasta",sep=""))
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

