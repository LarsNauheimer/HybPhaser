#########################################
### Preparation of the phasing script ###
#########################################

prep_phasing <- read.csv(csv_file_with_phasing_prep_info, header = T)
prep_phasing[is.na(prep_phasing)] <- ""
refseqs <- list.files(file.path(path_to_HybPhaser_results, "sequence_lists", reference_sequence_folder))
refseqs_fullpath <- list.files(path_to_reference_sequences, full.names = T)
reads <- list.files(path_to_read_files_phasing, full.names = T)

dir.create(folder_for_phased_reads, recursive = TRUE, showWarnings = FALSE)
dir.create(folder_for_phasing_stats, recursive = TRUE, showWarnings = FALSE)

if(!is.numeric(no_of_threads_phasing)){no_of_threads_phasing=1} 

# prepare reference part of the command

ref_command = vector()

for(i in 1:length(rownames(prep_phasing))){
  
  nrefs <- (length(which(prep_phasing[i,] !=""))-1)/2
  ref_co = vector()

  for(r in 1:nrefs){
    ref_file <- grep(prep_phasing[i,r*2],refseqs_fullpath, value=T)
    ref_abb <-prep_phasing[i,r*2+1] 
    ref_co[r] <- (paste("ref_",ref_abb,"=",ref_file, sep=""))
  }
ref_command[i] <- paste(ref_co, collapse = " ")
}  
  

# preparing the whole command line

phasing_command <- vector()

if(read_type_4phasing == "paired-end"){
  
  for(i in 1:length(rownames(prep_phasing))){
    
    read_file_1 <- reads[grep(paste(prep_phasing[i,1],ID_read_pair1,sep=""),reads)]
    read_file_2 <- reads[grep(paste(prep_phasing[i,1],ID_read_pair2,sep=""),reads)]
    
    phasing_command[i] <- paste(path_to_bbmap_executables,"bbsplit.sh ambiguous=all ambiguous2=all threads=",no_of_threads_phasing," ",ref_command[i], " in=",read_file_1," in2=",read_file_2," basename=",folder_for_phased_reads,prep_phasing[i,1],"_to_%.fastq", " refstats=",folder_for_phasing_stats,prep_phasing[i,1],"_phasing-stats.txt", sep="")
    
  }  

} else if (read_type_4phasing == "single-end") {

  for(i in 1:length(rownames(prep_phasing))){
    
    read_file <- reads[grep(paste(prep_phasing[i,1]),reads)]
    
    if(length(read_file)>1) { print("ERROR: multiple single end reads are selected!")}
    
    phasing_command[i] <- paste(path_to_bbmap_executables,"bbsplit.sh ambiguous=all ambiguous2=all threads=",no_of_threads_phasing," ",ref_command[i], " in=",read_file, " basename=",folder_for_phased_reads,prep_phasing[i,1],"_to_%.fastq", " refstats=",folder_for_phasing_stats,prep_phasing[i,1],"_phasing-stats.txt", sep="")
    
  }    

} else { 
  
  print("ERROR READ TYPE NOT SELECTED")
  
}






#### preparing phasing bash script
         
phasing_script_file <- file.path(path_to_phasing_folder,"run_bbsplit4phasing.sh")           
write("#!/bin/bash",file=phasing_script_file, append=F)
write(phasing_command, file=phasing_script_file, append=T)           
system(command = paste("chmod +x",phasing_script_file))
           
# run in R if selected
if(run_bash_script_in_R=="yes"){system(command = phasing_script_file)}

           
           
        
