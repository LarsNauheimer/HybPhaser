####################################
### Generation of sequence lists ###
####################################

# load config
if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)
library(seqinr)
library(stringr)

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_particular_dataset_optimization_subset)
output_sequences <- file.path(path_to_output_folder,paste("03_sequence_lists", name_for_particular_dataset_optimization_subset, sep=""))

if(intronerated_contig=="yes"){
  intronerated_name <- "intronerated" 
  intronerated_underscore <- "_"
} else {
  intronerated_name <- ""
  intronerated_underscore <- ""
}

targets <- read.fasta(fasta_file_with_targets, as.string=TRUE, set.attributes = FALSE)
targets_name <- unique(gsub(".*-","",labels(targets)))
samples <- readLines(path_to_namelist)


outsamples_missing <- readRDS(file=file.path(output_Robjects,"outsamples_missing.Rds"))
outloci_missing <- readRDS(file=file.path(output_Robjects,"outloci_missing.Rds"))
outloci_para_all <- readRDS(file=file.path(output_Robjects,"outloci_para_all.Rds"))
outloci_para_each <- readRDS(file=file.path(output_Robjects,"outloci_para_each.Rds"))
tab_snps_cl2b <- readRDS(file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))


#############################################


folder4seq_consensus_loci_raw <- file.path(output_sequences,"loci_raw_consensus")
folder4seq_consensus_loci_clean <- file.path(output_sequences,"loci_clean_consensus")
folder4seq_consensus_samples_raw <- file.path(output_sequences,"samples_raw_consensus")
folder4seq_consensus_samples_clean <- file.path(output_sequences,"samples_clean_consensus")
folder4seq_contig_loci_raw <- file.path(output_sequences,"loci_raw_contigs")
folder4seq_contig_loci_clean <- file.path(output_sequences,"loci_clean_contigs")
folder4seq_contig_samples_raw <- file.path(output_sequences,"samples_raw_contigs")
folder4seq_contig_samples_clean <- file.path(output_sequences,"samples_clean_contigs")

unlink(c(folder4seq_consensus_loci_raw,folder4seq_consensus_loci_clean, folder4seq_consensus_samples_raw, folder4seq_consensus_samples_clean, folder4seq_contig_loci_raw, folder4seq_contig_loci_clean, folder4seq_contig_samples_raw, folder4seq_contig_samples_clean),recursive = T) # delete directory, if it existed in order to prevent errors

dir.create(output_sequences, showWarnings = F)
dir.create(folder4seq_consensus_loci_raw, showWarnings = F)
dir.create(folder4seq_consensus_loci_clean, showWarnings = F)
dir.create(folder4seq_consensus_samples_raw, showWarnings = F)
dir.create(folder4seq_consensus_samples_clean, showWarnings = F)
dir.create(folder4seq_contig_loci_raw, showWarnings = F)
dir.create(folder4seq_contig_loci_clean, showWarnings = F)
dir.create(folder4seq_contig_samples_raw, showWarnings = F)
dir.create(folder4seq_contig_samples_clean, showWarnings = F)




#########################################################################
### concatenate consensus files across all samples to lists per locus ###
#########################################################################
# this extracts sequences from all subfolders in the HybPiper folder and collates them into one file per locus

# HybPhaser
if(intronerated_contig=="yes"){
  for(locus in targets_name){
    command_cat_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/*/intronerated_consensus/",paste(locus,"_intronerated.fasta",sep="")),">",file.path(folder4seq_consensus_loci_raw,paste(locus,"_intronerated_consensus_raw.fasta",sep="")))
    system(command_cat_consensus, ignore.stderr = TRUE)
    command_cat_contig <- paste("cat",file.path(path_to_output_folder,"01_data/*/intronerated_contigs/",paste(locus,"_intronerated.fasta",sep="")),">",file.path(folder4seq_contig_loci_raw,paste(locus,"_intronerated_contig_raw.fasta",sep="")))
    system(command_cat_contig, ignore.stderr = TRUE)
    command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci_raw,paste(locus,"_intronerated_consensus_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames_consensus) 
    command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci_raw,paste(locus,"_intronerated_contig_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames_contig)
  }
} else {
  for(locus in targets_name){
    command_cat_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/*/consensus/",paste(locus,".fasta",sep="")),">",file.path(folder4seq_consensus_loci_raw,paste(locus,"_consensus_raw.fasta",sep="")))
    command_cat_contig    <- paste("cat",file.path(path_to_output_folder,"01_data/*/contigs/",paste(locus,".fasta",sep="")),">",file.path(folder4seq_contig_loci_raw,paste(locus,"_contig_raw.fasta",sep="")))
    system(command_cat_consensus, ignore.stderr = TRUE)
    system(command_cat_contig, ignore.stderr = TRUE)
    command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci_raw,paste(locus,"_consensus_raw.fasta",sep="")), sep=""))
    command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci_raw,paste(locus,"_contig_raw.fasta",sep="")), sep=""))
    system(command_remove_locus_in_seqnames_consensus)
    system(command_remove_locus_in_seqnames_contig)
  } 
}  

# changing interleaved HybPiper files to non-interleaved fasta files
for(file in list.files(folder4seq_contig_loci_raw, full.names = T)){
  file <- list.files(folder4seq_contig_loci_raw, full.names = T)[3]
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}


# check whether accessions are in the hybpiper folder but not in the samples list.
# if the sample list is smaller some sequences have to be removed from the loci lists

hybpiper_result_dirs <- list.dirs(file.path(path_to_output_folder,"01_data"), full.names = FALSE, recursive = FALSE)
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
file.rename(list.files(folder4seq_consensus_loci_clean, full.names = T), gsub("_raw.fasta",".fasta",list.files(folder4seq_consensus_loci_clean, full.names = T)))




system(paste("for f in ",folder4seq_consensus_loci_clean,"/*_raw.fasta; do mv $f ${f/_raw/}; done", sep=""))

# contig
file.remove(list.files(folder4seq_contig_loci_clean,full.names = T)) # clear folder
file.copy(from=list.files(folder4seq_contig_loci_raw,full.names = T), to = folder4seq_contig_loci_clean) # copy raw files
file.rename(list.files(folder4seq_contig_loci_clean, full.names = T), gsub("_raw.fasta",".fasta",list.files(folder4seq_contig_loci_clean, full.names = T)))


## remove loci (failed/missing data/paralogs for all) for all samples

loci_files_consensus_clean <- list.files(path = folder4seq_consensus_loci_clean, full.names = T )
loci_files_contig_clean <- list.files(path = folder4seq_contig_loci_clean, full.names = T )

loci_to_remove <- c(names(failed_loci), outloci_missing, outloci_para_all)

if(length(loci_to_remove)!=0){
  loci_files_to_remove_consensus <- loci_files_consensus_clean[grep(paste(c(loci_to_remove),collapse="|"),loci_files_consensus_clean)]
  loci_files_to_remove_contig <- loci_files_contig_clean[grep(paste(c(loci_to_remove),collapse="|"),loci_files_contig_clean)]
  file.remove(loci_files_to_remove_consensus)
  file.remove(loci_files_to_remove_contig)
}


# get vector of all samples that should be removed from every locus

samples_to_remove_4all <- vector()

if(length(failed_samples) > 0 ){
  samples_to_remove_4all <- names(failed_samples)
} 

if(length(outsamples_missing) > 0 ){
  samples_to_remove_4all <-  unique(c(samples_to_remove_4all, outsamples_missing))
} 


## remove samples to be removed from all and sequences in each locus file from paralogs for each sample


for(locus in rownames(tab_snps_cl2b)){
  
  if(length(grep(locus,outloci_para_each)) >0 ){
    samples_to_remove <- c(samples_to_remove_4all, names(outloci_para_each[grep(locus,outloci_para_each)]))
  } else {
    samples_to_remove <- samples_to_remove_4all
  }
  
  
  if(length(samples_to_remove) !=0 ){
    
    # consensus
    locus_consensus_raw <- readLines(file.path(folder4seq_consensus_loci_clean,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
    lines_with_samplename <- grep(paste(samples_to_remove,collapse="|"),locus_consensus_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus_raw[-lines_to_remove]
      conn <- file(file.path(folder4seq_consensus_loci_clean,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
      writeLines(locus_file_red, conn)
      close(conn)
    }
    
    # contig
    locus_contig_raw <- readLines(file.path(folder4seq_contig_loci_clean,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
    lines_with_samplename <- grep(paste(samples_to_remove,collapse="|"),locus_contig_raw)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig_raw[-lines_to_remove]
      conn <- file(file.path(folder4seq_contig_loci_clean,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
    
  }
} 



#########################################################################
### concatenate consensus files across all loci to lists per sample   ###
#########################################################################

samples <- readLines(path_to_namelist)

# remove failed samples from list
if(length(failed_samples) != 0){
  samples <- samples[-which(samples %in% names(failed_samples))]
}

# collect all sequences

if(intronerated_contig=="yes"){
  for(sample in samples){
    command_cat_loci_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/intronerated_consensus/*.fasta"),">",file.path(folder4seq_consensus_samples_raw,paste(sample,"_intronerated_consensus_raw.fasta",sep="")))
    system(command_cat_loci_consensus)
    command_cat_loci_contig <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/intronerated_contigs/*.fasta"),">",file.path(folder4seq_contig_samples_raw,paste(sample,"_intronerated_contig_raw.fasta",sep="")))
    system(command_cat_loci_contig)
  } 
} else {
  for(sample in samples){
    command_cat_loci_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/consensus/*.fasta"),">",file.path(folder4seq_consensus_samples_raw,paste(sample,"_consensus_raw.fasta",sep="")))
    system(command_cat_loci_consensus)
    command_cat_loci_contig <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/contigs/*.fasta"),">",file.path(folder4seq_contig_samples_raw,paste(sample,"_contig_raw.fasta",sep="")))
    system(command_cat_loci_contig)
  }
}  


# changing interleaved HybPiper files to non-interleaved fasta files
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
file.rename(list.files(folder4seq_consensus_samples_clean, full.names = T), gsub("_raw.fasta",".fasta",list.files(folder4seq_consensus_samples_clean, full.names = T)))

#contig
file.remove(list.files(folder4seq_contig_samples_clean,full.names = T)) # clear folder
file.copy(from=list.files(folder4seq_contig_samples_raw,full.names = T), to = folder4seq_contig_samples_clean)
file.rename(list.files(folder4seq_contig_samples_clean, full.names = T), gsub("_raw.fasta",".fasta",list.files(folder4seq_contig_samples_clean, full.names = T)))



## remove all outlier loci (empty and high allele divergence) (cleaning step 1 and 2)
sample_files_consensus_clean <- list.files(path = folder4seq_consensus_samples_clean, full.names = T )
sample_files_contig_clean <- list.files(path = folder4seq_contig_samples_clean, full.names = T )


if(length(outsamples_missing)==0){
  samples_files_to_remove_consensus=""
  samples_files_to_remove_contig=""
} else {
  samples_files_to_remove_consensus <- sample_files_consensus_clean[grep(paste(c(outsamples_missing),collapse="|"),sample_files_consensus_clean)]
  samples_files_to_remove_contig <- sample_files_contig_clean[grep(paste(c(outsamples_missing),collapse="|"),sample_files_contig_clean)]
  file.remove(samples_files_to_remove_consensus)
  file.remove(samples_files_to_remove_contig)
}




samples_in <- samples
if(length(outsamples_missing) != 0){
  samples_in <- samples_in[-which(samples %in% outsamples_missing)]
}


loci_to_remove_4all <- vector()

if(length(failed_loci) > 0){
  loci_to_remove_4all <- names(failed_loci)
}

if(length(outloci_missing) > 0){
  loci_to_remove_4all <- c(loci_to_remove_4all, outloci_missing)
}

if(length(outloci_para_all) > 0){
  loci_to_remove_4all <- c(loci_to_remove_4all, outloci_para_all)
}

# remove outlier loci per sample in sample lists
for(sample in samples_in){
  
  
  if(length(grep(sample,names(outloci_para_each))) > 0 ){
    loci_to_remove <- c(loci_to_remove_4all, names(outloci_para_each[[which(names(outloci_para_each)%in% sample)]]))
  } else {
    loci_to_remove <- loci_to_remove_4all
  }
  
  
  if(length(loci_to_remove)!=0){
    
    
    #consensus
    consensus_file2clean <- file.path(folder4seq_consensus_samples_clean, paste(sample, intronerated_underscore, intronerated_name,"_consensus.fasta",sep=""))
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
    contig_file2clean <- file.path(folder4seq_contig_samples_clean, paste(sample,intronerated_underscore, intronerated_name,"_contig.fasta",sep=""))
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

