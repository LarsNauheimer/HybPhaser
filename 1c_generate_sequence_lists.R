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

if(name_for_dataset_optimization_subset != ""){
  folder_subset_add <- paste("_",name_for_dataset_optimization_subset, sep="")
} else {
  folder_subset_add <- ""
} 

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
output_sequences <- file.path(path_to_output_folder,paste("03_sequence_lists", folder_subset_add, sep=""))

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

tab_snps <- as.matrix(tab_snps_cl2b)
loci <- t(tab_snps)
failed_loci <- which(colSums(is.na(loci))==nrow(loci))
failed_samples <- which(colSums(is.na(tab_snps))==nrow(tab_snps))


#############################################

folder4seq_consensus_loci <- file.path(output_sequences,"loci_consensus")
folder4seq_consensus_samples <- file.path(output_sequences,"samples_consensus")
folder4seq_contig_loci <- file.path(output_sequences,"loci_contigs")
folder4seq_contig_samples <- file.path(output_sequences,"samples_contigs")


unlink(c(folder4seq_consensus_loci, folder4seq_consensus_samples,  folder4seq_contig_loci, folder4seq_contig_samples),recursive = T) # delete directory, if it existed in order to prevent errors

dir.create(output_sequences, showWarnings = F)
dir.create(folder4seq_consensus_loci, showWarnings = F)
dir.create(folder4seq_consensus_samples, showWarnings = F)
dir.create(folder4seq_contig_loci, showWarnings = F)
dir.create(folder4seq_contig_samples, showWarnings = F)




#########################################################################
### concatenate consensus files across all samples to lists per locus ###
#########################################################################
# this extracts sequences from all subfolders in the HybPiper folder and collates them into one file per locus


# if you are on Linux, then the bash commands cat and sed are used. If not, then R file operations are used, which are slower

if(Sys.info()['sysname']=="Linux"){
  if(intronerated_contig=="yes"){
    for(locus in targets_name){
      command_cat_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/*/intronerated_consensus/",paste(locus,"_intronerated.fasta",sep="")),">",file.path(folder4seq_consensus_loci,paste(locus,"_intronerated_consensus.fasta",sep="")))
      system(command_cat_consensus, ignore.stderr = TRUE)
      command_cat_contig <- paste("cat",file.path(path_to_output_folder,"01_data/*/intronerated_contigs/",paste(locus,"_intronerated.fasta",sep="")),">",file.path(folder4seq_contig_loci,paste(locus,"_intronerated_contig.fasta",sep="")))
      system(command_cat_contig, ignore.stderr = TRUE)
      command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci,paste(locus,"_intronerated_consensus.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames_consensus) 
      command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci,paste(locus,"_intronerated_contig.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames_contig)
    }
  } else {
    for(locus in targets_name){
      command_cat_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/*/consensus/",paste(locus,".fasta",sep="")),">",file.path(folder4seq_consensus_loci,paste(locus,"_consensus.fasta",sep="")))
      command_cat_contig    <- paste("cat",file.path(path_to_output_folder,"01_data/*/contigs/",paste(locus,".fasta",sep="")),">",file.path(folder4seq_contig_loci,paste(locus,"_contig.fasta",sep="")))
      system(command_cat_consensus, ignore.stderr = TRUE)
      system(command_cat_contig, ignore.stderr = TRUE)
      command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci,paste(locus,"_consensus.fasta",sep="")), sep=""))
      command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci,paste(locus,"_contig.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames_consensus)
      system(command_remove_locus_in_seqnames_contig)
    } 
  }  
} else {
  if(intronerated_contig=="yes"){
    for(locus in targets_name){
      #list all fasta files from that locus for all samples
      fasta_files <- list.files(path=file.path(path_to_output_folder,"01_data/"), pattern=paste(locus,"_intronerated.fasta",sep=""), recursive=TRUE, full.names = TRUE)
      #select consensus/contig files
      consensus_files <- grep("consensus",fasta_files,value = TRUE)
      contigs_files <- grep("contigs",fasta_files,value = TRUE)
      #define output files
      output_file_consensus <- file.path(folder4seq_consensus_loci,paste(locus,"_intronerated_consensus.fasta",sep=""))
      output_file_contigs <- file.path(folder4seq_contig_loci,paste(locus,"_intronerated_contigs.fasta",sep=""))
      #generate output files
      file.create(output_file_consensus, overwrite=TRUE)
      file.create(output_file_contigs, overwrite=TRUE)
      #append sample fastas to the empty output file
      file.append(output_file_consensus,consensus_files)
      file.append(output_file_contigs,contigs_files)
      #read lines of each file and remove the "-locus" of the sequence names
      lines_consensus <- gsub(paste("-",locus,sep=""),"",readLines(output_file_consensus))
      lines_contigs <- gsub(paste("-",locus,sep=""),"",readLines(output_file_contigs))
      #write lines into files
      write(lines_consensus, file=output_file_consensus)
      write(lines_contigs, file=output_file_contigs)
    }
  } else {
    for(locus in targets_name){
      #list all fasta files from that locus for all samples
      fasta_files <- list.files(path=file.path(path_to_output_folder,"01_data/"), pattern=paste(locus,".fasta",sep=""), recursive=TRUE, full.names = TRUE)
      #select consensus/contig files
      consensus_files <- grep("consensus",fasta_files,value = TRUE)
      contigs_files <- grep("contigs",fasta_files,value = TRUE)
      #define output files
      output_file_consensus <- file.path(folder4seq_consensus_loci,paste(locus,"_consensus.fasta",sep=""))
      output_file_contigs <- file.path(folder4seq_contig_loci,paste(locus,"_contigs.fasta",sep=""))
      #generate output files
      file.create(output_file_consensus, overwrite=TRUE)
      file.create(output_file_contigs, overwrite=TRUE)
      #append sample fastas to the empty output file
      file.append(output_file_consensus,consensus_files)
      file.append(output_file_contigs,contigs_files)
      #read lines of each file and remove the "-locus" of the sequence names
      lines_consensus <- gsub(paste("-",locus,sep=""),"",readLines(output_file_consensus))
      lines_contigs <- gsub(paste("-",locus,sep=""),"",readLines(output_file_contigs))
      #write lines into files
      write(lines_consensus, file=output_file_consensus)
      write(lines_contigs, file=output_file_contigs)
    } 
  }  
}


# changing interleaved HybPiper files to non-interleaved fasta files
for(file in list.files(folder4seq_contig_loci, full.names = T)){
  file <- list.files(folder4seq_contig_loci, full.names = T)[3]
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
  for(raw_consensus_file in list.files(folder4seq_consensus_loci, full.names = TRUE)){
    locus_consensus <- readLines(raw_consensus_file)
    lines_with_samplename <- which(gsub(">","",locus_consensus) %in% dirs_not_in_sample_list)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus[-lines_to_remove]
      conn <- file(raw_consensus_file)
      writeLines(locus_file_red, conn)
      close(conn)
    }
  }
  # contig
  for(raw_contig_file in list.files(folder4seq_contig_loci, full.names = TRUE)){
    locus_contig <- readLines(raw_contig_file)
    lines_with_samplename <- which(gsub(">","",locus_contig) %in% dirs_not_in_sample_list)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig[-lines_to_remove]
      conn <- file(raw_contig_file)
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
  }
}


### remove loci from dataset optimization (missing data and paralogs)
######################################################################

## remove loci (failed/missing data/paralogs for all) for all samples

loci_files_consensus <- list.files(path = folder4seq_consensus_loci, full.names = T )
loci_files_contig <- list.files(path = folder4seq_contig_loci, full.names = T )

loci_to_remove <- c(names(failed_loci), outloci_missing, outloci_para_all)

if(length(loci_to_remove)!=0){
	if(intronerated_contig=="no"){
		loci_files_to_remove_consensus <- loci_files_consensus[which(gsub(".*/(.*)_consensus.fasta","\\1",loci_files_consensus) %in% loci_to_remove)]
		loci_files_to_remove_contig <- loci_files_contig[which(gsub(".*/(.*)_contig.fasta","\\1",loci_files_contig) %in% loci_to_remove)]
	} else {
		loci_files_to_remove_consensus <- loci_files_consensus[which(gsub(".*/(.*)_intronerated_consensus.fasta","\\1",loci_files_consensus) %in% loci_to_remove)]
		loci_files_to_remove_contig <- loci_files_contig[which(gsub(".*/(.*)_intronerated_contig.fasta","\\1",loci_files_contig) %in% loci_to_remove)]
	}
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
  if(!(locus %in% names(failed_loci))){
  
    if(length(grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each)) >0 ){
      samples_to_remove <- c(samples_to_remove_4all, names(outloci_para_each[grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each)]))
    } else {
      samples_to_remove <- samples_to_remove_4all
    }
    
    
    if(length(samples_to_remove) !=0 ){
      
      # consensus
      locus_consensus <- readLines(file.path(folder4seq_consensus_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
      lines_with_samplename <- which(gsub(">","",locus_consensus) %in% samples_to_remove)
      if(length(lines_with_samplename) !=0){
        lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
        locus_file_red <- locus_consensus[-lines_to_remove]
        conn <- file(file.path(folder4seq_consensus_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
        writeLines(locus_file_red, conn)
        close(conn)
      }
      
      # contig
      locus_contig <- readLines(file.path(folder4seq_contig_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
      lines_with_samplename <- which(gsub(">","",locus_contig) %in% samples_to_remove)
      if(length(lines_with_samplename) !=0){
        lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
        locus_file_hp_red <- locus_contig[-lines_to_remove]
        conn <- file(file.path(folder4seq_contig_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
        writeLines(locus_file_hp_red, conn)
        close(conn)
      }
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
# if you are on Linux, then the bash commands cat and sed are used. If not, then R file operations are used, which are slower
if(Sys.info()['sysname']=="Linux"){
  if(intronerated_contig=="yes"){
    for(sample in samples){
      command_cat_loci_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/intronerated_consensus/*.fasta"),">",file.path(folder4seq_consensus_samples,paste(sample,"_intronerated_consensus.fasta",sep="")))
      system(command_cat_loci_consensus)
      command_cat_loci_contig <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/intronerated_contigs/*.fasta"),">",file.path(folder4seq_contig_samples,paste(sample,"_intronerated_contig.fasta",sep="")))
      system(command_cat_loci_contig)
    } 
  } else {
    for(sample in samples){
      command_cat_loci_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/consensus/*.fasta"),">",file.path(folder4seq_consensus_samples,paste(sample,"_consensus.fasta",sep="")))
      system(command_cat_loci_consensus)
      command_cat_loci_contig <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/contigs/*.fasta"),">",file.path(folder4seq_contig_samples,paste(sample,"_contig.fasta",sep="")))
      system(command_cat_loci_contig)
    }
  }  
} else {
  if(intronerated_contig=="yes"){
    for(sample in samples){
      #list all fasta files from the sample for all loci
      consensus_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/intronerated_consensus/"),pattern="*.fasta", full.names = TRUE)
      contigs_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/intronerated_contigs/"),pattern="*.fasta", full.names = TRUE)
      #define output files
      output_file_consensus_samples <- file.path(folder4seq_consensus_samples,paste(sample,"_intronerated_consensus.fasta",sep=""))          
      output_file_contigs_samples <- file.path(folder4seq_contig_samples,paste(sample,"_intronerated_contigs.fasta",sep=""))          
      #create output files
      file.create(output_file_consensus_samples, overwrite=TRUE)
      file.create(output_file_contigs_samples, overwrite=TRUE)
      #append fasta files to empty output file
      file.append(output_file_consensus_samples,consensus_files_samples)                
      file.append(output_file_contigs_samples,contigs_files_samples)        
    } 
  } else {
    for(sample in samples){
      
      #list all fasta files from the sample for all loci
      consensus_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/consensus/"),pattern="*.fasta", full.names = TRUE)
      contigs_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/contigs/"),pattern="*.fasta", full.names = TRUE)
      #define output files
      output_file_consensus_samples <- file.path(folder4seq_consensus_samples,paste(sample,"_consensus.fasta",sep=""))          
      output_file_contigs_samples <- file.path(folder4seq_contig_samples,paste(sample,"_contigs.fasta",sep=""))          
      #create output files
      file.create(output_file_consensus_samples, overwrite=TRUE)
      file.create(output_file_contigs_samples, overwrite=TRUE)
      #append fasta files to empty output file
      file.append(output_file_consensus_samples,consensus_files_samples)                
      file.append(output_file_contigs_samples,contigs_files_samples)        
    }
  }
}

# changing interleaved HybPiper files to non-interleaved fasta files
for(file in list.files(folder4seq_contig_samples, full.names = T)){
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}




# remove samples (missing data, paralogs for all)
###################################################

## copy all sequence lists to new folder before removing parts of it


## remove all outlier loci (empty and high allele divergence) (cleaning step 1 and 2)
sample_files_consensus <- list.files(path = folder4seq_consensus_samples, full.names = T )
sample_files_contig <- list.files(path = folder4seq_contig_samples, full.names = T )


if(length(outsamples_missing)==0){
  samples_files_to_remove_consensus=""
  samples_files_to_remove_contig=""
} else {
  samples_files_to_remove_consensus <- sample_files_consensus[which(gsub(".*/(.*)_consensus.fasta","\\1",sample_files_consensus) %in% outsamples_missing)]
  samples_files_to_remove_contig <- sample_files_contig[which(gsub(".*/(.*)_contigs.fasta","\\1",sample_files_contig) %in% outsamples_missing)]

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
  
  
  if(length(which(names(outloci_para_each)%in% sample)) > 0 ){
    loci_to_remove <- c(loci_to_remove_4all, names(outloci_para_each[[which(names(outloci_para_each)%in% sample)]]))
  } else {
    loci_to_remove <- loci_to_remove_4all
  }
  
  
  
  if(length(loci_to_remove)!=0){
    
    
    #consensus
    consensus_file2clean <- file.path(folder4seq_consensus_samples, paste(sample, intronerated_underscore, intronerated_name,"_consensus.fasta",sep=""))
    samples_consensus <- readLines(consensus_file2clean)
    lines_with_lociname <- which(gsub(">.*-","",samples_consensus) %in% loci_to_remove)
    if(length(lines_with_lociname) !=0){
      lines_to_remove <- c(lines_with_lociname,lines_with_lociname+1)
      sample_file_consensus_red <- samples_consensus[-lines_to_remove]
      conn <- file(consensus_file2clean)
      writeLines(sample_file_consensus_red, conn)
      close(conn)
    }
    
    #contig
    contig_file2clean <- file.path(folder4seq_contig_samples, paste(sample,intronerated_underscore, intronerated_name,"_contig.fasta",sep=""))
    samples_contig <- readLines(contig_file2clean)
    lines_with_lociname <- which(gsub(">.*-","",samples_contig) %in% loci_to_remove)
    if(length(lines_with_lociname) !=0){
      lines_to_remove <- c(lines_with_lociname,lines_with_lociname+1)
      sample_file_contig_red <- samples_contig[-lines_to_remove]
      conn <- file(contig_file2clean)
      writeLines(sample_file_contig_red, conn)
      close(conn)
    }
  }
} 
