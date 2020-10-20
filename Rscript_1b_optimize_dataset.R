################################################################################
### Applying cleaning steps to 1) reduce missing data and 2) remove paralogs ###
################################################################################

# retrieving data from previous steps 
tab_snps <- readRDS(file=file.path(output_Robjects,"Table_SNPs_raw.Rds"))
tab_length <- readRDS(file=file.path(output_Robjects,"Table_consensus_length.Rds"))

rownames(tab_snps) <-  tab_snps[,1]
tab_snps <- as.matrix(tab_snps[,-c(1)])
loci <- t(tab_snps)

dir.create(output_cleaning, showWarnings = F)

#####################################
### cleaning step 1: missing data ### 
#####################################

# 1a) Failed loci and samples
#############################

# checking for and removing loci without any sequences
failed_loci <- which(colSums(is.na(loci))==nrow(loci))
loci <- loci[,colSums(is.na(loci))<nrow(loci)]
tab_snps <- tab_snps[rowSums(is.na(tab_snps))<ncol(tab_snps),]


# checking for and removing samples without any sequences
failed_samples <- which(colSums(is.na(tab_snps))==nrow(tab_snps))
tab_snps <- tab_snps[,colSums(is.na(tab_snps))<nrow(tab_snps)]

# output
summary_file=file.path(output_cleaning,"Summary_missing_data.txt")
cat(file=summary_file,"Summary of cleaning step 1",paste("There were: ",length(failed_loci), "failed loci:",sep=""),names(failed_loci),paste("There were ",length(failed_samples), " failed samples:",sep=""),names(failed_samples),"", sep="\n")


### 1b) Loci with missing data
#################################

nloci <- length(colnames(loci))
nsamples <- length(rownames(loci))

loci_missing <- vector()
for(i in 1:nloci){
  loci_missing[i] <- length(which(is.na(loci[,i]))) 
}
names(loci_missing) <- colnames(loci)

# generate PDF and PNG files
pdf(file=file.path(output_cleaning,"1b_Missing_data_per_locus.pdf"), width = 8, height=6)
  par(mfrow=c(1,2))
  boxplot(loci_missing, main="Missing data per locus", sub=paste("loci=",nloci, " samples=",nsamples))
  hist(loci_missing, breaks=50, main="")
dev.off()

png(file=file.path(output_cleaning,"1b_Missing_data_per_locus.png"), width = 800, height=600)
  par(mfrow=c(1,2))
  boxplot(loci_missing, main="Missing data per locus", sub=paste("loci=",nloci, " samples=",nsamples))
  hist(loci_missing, breaks=50, main="")
dev.off()

par(mfrow=c(1,1))


if (length(threshold_exclude_loci_x_missing_samples) == 0 || threshold_exclude_loci_x_missing_samples == "none"){
  outloci_missing <-  vector()
} else if (threshold_exclude_loci_x_missing_samples == "1.5*IQR"){
  outloci_missing <- boxplot(loci_missing, plot=F) $out
  outloci_missing <- outloci_missing[which(outloci_missing > median(loci_missing))]  # only use high outliers
} else if (threshold_exclude_loci_x_missing_samples == "proportion"){
  outloci_missing <- loci_missing[which(loci_missing >= proportion_of_missing_samples_per_locus * nsamples)] 
} else {
  outloci_missing <- loci_missing[which(loci_missing >= threshold_exclude_loci_x_missing_samples)] 
}


#output
cat(file=summary_file, append = T, "\nCleaning step 1b: missing data (samples) of loci\n")
cat(file=summary_file, append = T, "The average number of missing samples for all loci is:", round(mean(loci_missing),2), "(of",nsamples, "samples)\n")
cat(file=summary_file, append = T, "The median number of missing samples for all loci is:", median(loci_missing), "(of",nsamples, "samples)\n")
if(threshold_exclude_loci_x_missing_samples == "proportion"){
  cat(file=summary_file, append = T, "The chosen threshold for removing loci with missing data is set as the proportion (",proportion_of_missing_samples_per_locus," of", nsamples," samples).\n")
} else {
  cat(file=summary_file, append = T, "The chosen threshold for removing loci with missing data is:",threshold_exclude_loci_x_missing_samples,"\n")
}
cat(file=summary_file, append = T, length(outloci_missing),"loci will be removed:\n", paste(names(outloci_missing),outloci_missing ,"\n"))

#write.csv( data.frame( "# of empty samples"=outloci_missing), file.path(output_cleaning, "1b_loci_removed_for_missing_data.csv"))



### 1c) Samples with missing data
####################################

samples_missing <- 0
for(i in 1:length(colnames(tab_snps))){
  samples_missing[i] <- length(which(is.na(tab_snps[,i]))) 
}
names(samples_missing) <- colnames(tab_snps)

# graphs to check: boxplot showing the distribution and density plot  
# generate PDF and PNG files 
pdf(file=file.path(output_cleaning,"1c_Missing_data_per_sample.pdf"), width = 8, height=6)
 par(mfrow=c(1,2))
 boxplot(samples_missing, main="Missing data per sample", sub=paste("loci=",nloci, " samples=",nsamples))
 hist(samples_missing, breaks=50, main="")
dev.off()

png(file=file.path(output_cleaning,"1c_Missing_data_per_sample.png"), width = 800, height=600)
 par(mfrow=c(1,2))
 boxplot(samples_missing, main="Missing data per sample", sub=paste("loci=",nloci, " samples=",nsamples))
 hist(samples_missing, breaks=50, main="")
dev.off()

par(mfrow=c(1,1))


if (length(threshold_exclude_samples_x_missing_loci) == 0 || threshold_exclude_samples_x_missing_loci == "none"){
  outsamples_missing <-  vector()
} else if (threshold_exclude_samples_x_missing_loci == "1.5*IQR"){
  outsamples_missing <- boxplot(samples_missing, plot=F) $out
  outsamples_missing <- outsamples_missing[which(outsamples_missing > median(samples_missing))]  # only use high outliers
} else if (threshold_exclude_samples_x_missing_loci == "proportion"){
  outsamples_missing <- samples_missing[which(samples_missing >= proportion_of_missing_loci_per_sample * nloci)] 
} else {
  outsamples_missing <- samples_missing[which(samples_missing >= threshold_exclude_samples_x_missing_loci)] 
}


# ouput to file
cat(file=summary_file, append = T, "\nCleaning step 1c: missing data (loci) of samples\n")
cat(file=summary_file, append = T, "The average number of missing loci for all samples is:", round(mean(samples_missing),2), "(of",nloci, "loci)\n")
cat(file=summary_file, append = T, "The median number of missing loci for all samples is:", median(samples_missing), "(of",nloci, "loci)\n")
if(threshold_exclude_samples_x_missing_loci == "proportion"){
  cat(file=summary_file, append = T, "The chosen threshold for removing samples with missing data is set as the proportion (",proportion_of_missing_loci_per_sample," of", nloci," loci).\n")
} else {
  cat(file=summary_file, append = T, "The chosen threshold for removing samples with missing data is:",threshold_exclude_samples_x_missing_loci,"\n")
}
cat(file=summary_file, append = T, length(outsamples_missing),"samples will be removed:\n", paste(names(outsamples_missing),outsamples_missing ,"\n"))



### 1d) Removing samples with poor sequence recovery (in bp as proportion of the target sequence lengths)
#########################################################################################################

tab_length <- readRDS(file=file.path(output_Robjects,"Table_consensus_length.Rds"))

targets_length_all <- lengths(read.FASTA(fasta_file_with_targets))
gene_names <- unique(gsub(".*-","",names(targets_length_all)))
max_target_length <- vector()

for(i in 1:length(gene_names)){
  max_target_length[i] <- max(targets_length_all[grep(gene_names[i],names(targets_length_all))])
}

comb_target_length <- sum(max_target_length)

comb_seq_length <- colSums(tab_length[,-1], na.rm = T)
prop_of_target <- comb_seq_length/comb_target_length

# graphs to check: boxplot showing the distribution and density plot  
# generate PDF and PNG files 
pdf(file=file.path(output_cleaning,"1d_Sequence_length_as_proportion_of_combined_targets_length.pdf"), width = 8, height=6)
  par(mfrow=c(1,2))
  boxplot(prop_of_target, main="Sequence length as prop.\nof comb. targets length", sub=paste("combined target length:\n",comb_target_length,"bp" ))
  hist(prop_of_target, breaks=50, main="")
dev.off()

png(file=file.path(output_cleaning,"1d_Total_recovered_sequence_length_as_proportion_of_targets_length.png"), width = 800, height=600)
  par(mfrow=c(1,2))
  boxplot(prop_of_target, main="Sequence length as prop.\nof comb. targets length", sub=paste("combined target length:\n",comb_target_length,"bp" ))
  hist(prop_of_target, breaks=50, main="")
dev.off()

par(mfrow=c(1,1))


if (length(threshold_exclude_samples_x_proportion_of_target_sequence) == 0 || threshold_exclude_samples_x_proportion_of_target_sequence == "none"){
  outsamples_recovered_seq_length <-  vector()
} else {
  outsamples_recovered_seq_length <- prop_of_target[which(prop_of_target <= threshold_exclude_samples_x_proportion_of_target_sequence)] 
}


# ouput to file
cat(file=summary_file, append = T, "\nCleaning step 1d: Combined sequence length recovered for each sample\n")
cat(file=summary_file, append = T, "The total length of combined target sequences is: ", comb_target_length," bp.")
cat(file=summary_file, append = T, "The average proportion of targets length recovered is:", round(mean(prop_of_target),2),".")
cat(file=summary_file, append = T, "The median number of missing loci for all samples is:",  round(median(prop_of_target),2),".")
cat(file=summary_file, append = T, "The chosen threshold for proportion of combined target length to exclude low recovery samples: ",threshold_exclude_samples_x_missing_loci,".\n")
cat(file=summary_file, append = T, length(outsamples_recovered_seq_length),"samples will be removed:\n", paste(names(outsamples_recovered_seq_length),outsamples_recovered_seq_length ,".\n"))





#########################################################
# removing bad loci and samples from the table 

tab_snps_cl1 <- tab_snps
if(length(outloci_missing) != 0){tab_snps_cl1 <- tab_snps_cl1[-which(rownames(tab_snps) %in% names(outloci_missing)),]}
if(length(outsamples_missing) != 0){tab_snps_cl1 <- tab_snps_cl1[,-which(colnames(tab_snps) %in% unique(names(c(outsamples_missing,outsamples_recovered_seq_length))))]}
loci_cl1 <- t(tab_snps_cl1)



################################################################################
### cleaning step 2, removing paralogs for a) all samples and b) each sample ### 
################################################################################


### 2a) Paralogs across multiple samples (removing outlier loci for all samples)
################################################################################

# (Loci that have a high proportion of SNPS across all samples are likely paralogs or contain other sources of contamination)

loci_cl1_colmeans <- colMeans(as.matrix(loci_cl1), na.rm = T)
nloci_cl1 <- length(colnames(loci_cl1))
nsamples_cl1 <- length(colnames(tab_snps_cl1))

loci_cl1_colmeans_mean <- round(mean(loci_cl1_colmeans),4)
loci_cl1_colmeans_median <- round(median(loci_cl1_colmeans),4)

# Distribution of loci and their mean proportion of SNPs
pdf(file=file.path(output_cleaning,"2a_Overview_SNPs_per_locus.pdf"), width = 8, height = 6)
  par(mfrow=c(1,2))
  boxplot(loci_cl1_colmeans, main="Mean % SNPs of loci", sub=paste("loci=",nloci_cl1, " samples=",nsamples_cl1))
  hist(loci_cl1_colmeans, breaks = 50, main='')
dev.off()
  
png(file=file.path(output_cleaning,"2a_Overview_SNPs_per_locus.png"), width = 800, height = 600)
  par(mfrow=c(1,2))
  boxplot(loci_cl1_colmeans, main="Mean % SNPs of loci", sub=paste("loci=",nloci_cl1, " samples=",nsamples_cl1))
  hist(loci_cl1_colmeans, breaks = 50, main='')
dev.off()

par(mfrow=c(1,1))



outloci_para_all <- boxplot(loci_cl1_colmeans, plot=F) $out
outloci_para_all <- outloci_para_all[which(outloci_para_all > median(loci_cl1_colmeans))]


if (length(threshold_exclude_loci_paralogs_4all_propSNP) == 0 || threshold_exclude_loci_paralogs_4all_propSNP == "none"){
  outloci_para_all <-   vector()
} else if (threshold_exclude_loci_paralogs_4all_propSNP == "1.5*IQR"){
  outloci_para_all <- boxplot(loci_cl1_colmeans) $out
  outloci_para_all <- outloci_para_all[which(outloci_para_all > median(loci_cl1_colmeans))]
} else {
  outloci_para_all <- loci_cl1_colmeans[which(loci_cl1_colmeans >= threshold_exclude_loci_paralogs_4all_propSNP)]
}


colour_outparaall <- rep("white",nloci_cl1)
colour_outparaall[which(colnames(loci_cl1[,order(loci_cl1_colmeans)]) %in% names(outloci_para_all))] <- "red"
loci_cl1_order_means <- loci_cl1[,order(loci_cl1_colmeans)]

pdf(file=file.path(output_cleaning,"2a_Boxplots_SNPs_per_locus.pdf"), width = 15, h=40)
  boxplot(loci_cl1_order_means, horizontal=T, col=colour_outparaall, las=2)
dev.off()

png(file=file.path(output_cleaning,"2a_Boxplots_SNPs_per_locus.png"), width = 1500, h=4000)
  boxplot(loci_cl1_order_means, horizontal=T, col=colour_outparaall, las=1)
dev.off()



### removing marked loci from table
if(length(outloci_para_all)==0) {tab_snps_cl2a <- tab_snps_cl1
} else { tab_snps_cl2a <- tab_snps_cl1[-which(rownames(tab_snps_cl1) %in%  names(outloci_para_all)),]}



### 2b) Paralogs for each sample (removing outlier loci for each sample)
##########################################################################


tab_snps_cl2b <- tab_snps_cl2a

#Check outlier per sample before cleaning3

boxplot(as.data.frame(tab_snps_cl2a[,order(colMeans(as.matrix(tab_snps_cl2a), na.rm = T))]), horizontal=T, las=1)

# generate graphs
pdf(file=file.path(output_cleaning,"2b_Boxplots_SNPs_per_sample_before_cleaning.pdf"), width = 15, h=30)
  boxplot(as.data.frame(tab_snps_cl2a[,order(colMeans(as.matrix(tab_snps_cl2a), na.rm = T))]), horizontal=T, las=1)
dev.off()

png(file=file.path(output_cleaning,"2b_Boxplots_SNPs_per_sample_before_cleaning.png"), width = 1500, h=3000)
  boxplot(as.data.frame(tab_snps_cl2a[,order(colMeans(as.matrix(tab_snps_cl2a), na.rm = T))]), horizontal=T, las=1)
dev.off()


if(!exists("remove_paralogs_4each")) {remove_paralogs_4each <- "no"}
  
if(remove_paralogs_4each == "yes" ){

  outloci_para_each <- list()
  outloci_para_each <- sapply(colnames(tab_snps_cl2a),function(x) NULL)
  
  for(i in 1:length(colnames(tab_snps_cl2a))){
    #i=15
    outlier_i <- boxplot(replace(tab_snps_cl2a[,i],tab_snps_cl2a[,i]==0, NA),plot=F) $out  # here we only regard heterozygous loci
    outlier_loci_i <- names(outlier_i)
    outloci_para_each[i] <- list(outlier_loci_i)
    #cat(file = outlier_gene_file,"Sample:",colnames(tab_snps_cl1)[i],"has", length(outlier_loci_i), "outlier genes (",outlier_loci_i,")\n", append=T)
    
    tab_snps_cl2b[which(rownames(tab_snps_cl2a) %in% outlier_loci_i),i] <- NA
  }
  

  
  # generate graphic files
  pdf(file=file.path(output_cleaning,"2b_Boxplots_SNPs_per_sample_clean.pdf"), width = 15, h=30)
    boxplot(as.data.frame(tab_snps_cl2b[,order(colMeans(as.matrix(tab_snps_cl2b), na.rm = T))]), horizontal=T, las=1)
  dev.off()
  
  png(file=file.path(output_cleaning,"2b_Boxplots_SNPs_per_sample_clean.png"), width = 1500, h=3000)
    boxplot(as.data.frame(tab_snps_cl2b[,order(colMeans(as.matrix(tab_snps_cl2b), na.rm = T))]), horizontal=T, las=1)
  dev.off()
  
} else {
  
  outloci_para_each <- list()
  outloci_para_each <- sapply(colnames(tab_snps_cl2a),function(x) NULL)
  tab_snps_cl2b <- tab_snps_cl2a
  
}

tab_length_cl2b <- tab_length[which(rownames(tab_length) %in% rownames(tab_snps_cl2b)),which(colnames(tab_length) %in% colnames(tab_snps_cl2b))]


#####################################
### output text files and r-data  ###
#####################################


# write text file with all loci that were outliers across all samples

cl2b_file <- file.path(output_cleaning,"Summary_Paralogs.txt")
cat(file=cl2b_file,"Summary of revoval of putative paralog loci")

cat(file=cl2b_file,"Paralogs removed for all samples:\n", append = T)
cat(file=cl2b_file, paste(names(outloci_para_all),"\n"), append = T)
cat(file=cl2b_file,  "\n\n", append = T)

if(length(outloci_para_each) > 0){
  cat(file=cl2b_file, "Paralogs removed for each sample:\n", append=T) 
  for(i in 1:length(names(outloci_para_each))){
    cat(file=cl2b_file, names(outloci_para_each)[i],"\t", append=T)  
    cat(file=cl2b_file, length(outloci_para_each[[i]]),outloci_para_each[[i]], sep="\t", append=T)  
    cat(file=cl2b_file, "\n", append=T)  
  }
} else{
  cat(file=cl2b_file, "The step for removing paralogs for each samples was skipped.\n", append=T) 
}

write.csv(tab_snps_cl2b,file=file.path(output_cleaning,"Table_SNPs_clean.csv"))
write.csv(tab_length_cl2b,file=file.path(output_cleaning,"Table_consensus_length_clean.csv"))

### save Data as R objects
output_Robjects <-file.path(path_for_HybPhaser_output,"R_objects")

saveRDS(tab_snps_cl2b,file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))
saveRDS(tab_length_cl2b,file=file.path(output_Robjects,"Table_consensus_length_cleaned.Rds"))

saveRDS(outloci_missing,file=file.path(output_Robjects,"outloci_missing.Rds"))
saveRDS(outsamples_missing,file=file.path(output_Robjects,"outsamples_missing.Rds"))
saveRDS(outsamples_recovered_seq_length,file=file.path(output_Robjects,"outsamples_recovered_seq_length.Rds"))
saveRDS(outloci_para_all,file=file.path(output_Robjects,"outloci_para_all.Rds"))
saveRDS(outloci_para_each,file=file.path(output_Robjects,"outloci_para_each.Rds"))

