###########################################################
### Preparation of BBSplit script for Clade association ###
###########################################################

# load config
if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)


if(read_type_cladeassociation == "paired-end"){
  
  read_files <- list.files(path_to_read_files_cladeassociation, full.names = F, pattern = paste("*",ID_read_pair1,".*", sep="") , include.dirs = F)  

} else if(read_type_cladeassociation == "single-end") {

  read_files <- list.files(path_to_read_files_cladeassociation, full.names = F, pattern = "*.*", include.dirs = F)  
  ID_read_pair1 <- ""
   
} else { 
  
  print("ERROR READ TYPE NOT SELECTED")

}

if(file_with_samples_included == "" | file_with_samples_included == "none" ){
   } else {
  samples_in <- readLines(file_with_samples_included)
  samples_in <- gsub("\\s*$","",samples_in)
  read_files_noID <- gsub(ID_read_pair1,"",read_files)
  read_files_noend <- gsub("[.]fast.*","",read_files_noID)
  read_files <- read_files[match(samples_in, read_files_noend)]
}


read_files <- na.exclude(read_files)

folder_bbsplit_stats <- file.path(path_to_clade_association_folder,"bbsplit_stats/")
dir.create(folder_bbsplit_stats, showWarnings = FALSE)

sample_sequences <- list.files(path_to_reference_sequences)

ref_samples <- read.csv(csv_file_with_clade_reference_names, header = T)
colnames(ref_samples)[1] <- "samples"
colnames(ref_samples)[2] <- "abb"

ref_samples_files <- sample_sequences[match(ref_samples$samples,gsub("(_intronerated|_consensus|_contig).*","",sample_sequences))]

if(length(ref_samples$abb)>0){
  ref_command <- paste(paste("ref_",ref_samples[,2],"=", path_to_reference_sequences,"/",ref_samples_files, sep=""), collapse=" ")
} else {
  ref_command <- paste("ref=",paste(path_to_reference_sequences,ref_samples_files, collapse = ",",sep=""),sep="")
}


# setting no of threads
if(no_of_threads_clade_association == 0 ||  no_of_threads_clade_association == "auto") {
  threadtext <- ""
} else {
  threadtext <- paste(" threads=",no_of_threads_clade_association,sep="")
}


# setting Java memory usage 

if(java_memory_usage_clade_association != ""){
  caXmx <- paste(" -Xmx",java_memory_usage_clade_association, sep="")
} else {
  caXmx <- ""
}

if(path_to_bbmap == ""){
  bbsplit_sh <- "bbsplit.sh"
} else {
  bbsplit_sh <- file.path(path_to_bbmap,"bbsplit.sh")
}



write("#!/bin/bash",file=file.path(path_to_clade_association_folder,"run_bbsplit4clade_association.sh"), append=F)
  
for(i in 1:length(read_files)){
    
    if(read_type_cladeassociation == "paired-end"){
      
      bbsplit_command <- paste(bbsplit_sh," ",ref_command, " in=",file.path(path_to_read_files_cladeassociation,read_files[i]), " in2=",sub(ID_read_pair1,ID_read_pair2,file.path(path_to_read_files_cladeassociation,read_files[i])),threadtext," ambiguous=random ambiguous2=all refstats=",folder_bbsplit_stats, gsub(ID_read_pair1,"",read_files[i]),"_bbsplit-stats.txt", caXmx, sep="")
      
    } else {
  
      bbsplit_command <- paste(bbsplit_sh," ",ref_command, " in=",file.path(path_to_read_files_cladeassociation,read_files[i]),threadtext," ambiguous=random ambiguous2=all refstats=",folder_bbsplit_stats,read_files[i],"_bbsplit-stats.txt", caXmx, sep="")
      
    }
  
    write(bbsplit_command,file=file.path(path_to_clade_association_folder,"run_bbsplit4clade_association.sh"), append=T)
  
}
  
system(command = paste("chmod +x",file.path(path_to_clade_association_folder,"run_bbsplit4clade_association.sh")))

if(run_clade_association_mapping_in_R=="yes"){
  system(command = file.path(path_to_clade_association_folder,"run_bbsplit4clade_association.sh"))
}
