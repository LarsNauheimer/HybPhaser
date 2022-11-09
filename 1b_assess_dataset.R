#####################################################################################
### Generating files for assessment of missing data,. paralogs and heterozygosity ###
#####################################################################################

# load config
if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)

# generate folders
if(name_for_dataset_optimization_subset != ""){
  folder_subset_add <- paste("_",name_for_dataset_optimization_subset, sep="")
}else {
  folder_subset_add <- ""
} 

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)

output_assess <- file.path(path_to_output_folder,paste("02_assessment", folder_subset_add, sep=""))

dir.create(output_assess, showWarnings = F)
dir.create(output_Robjects, showWarnings = F)

# check if SNPs count for subset has been done, if not use SNPs count of first run subset with name "")
if(file.exists(file=file.path(output_Robjects,"Table_SNPs.Rds"))){
  tab_snps <- readRDS(file=file.path(output_Robjects,"Table_SNPs.Rds"))
  tab_length <- readRDS(file=file.path(output_Robjects,"Table_consensus_length.Rds"))
} else {
  tab_snps <- readRDS(file=file.path(path_to_output_folder,"00_R_objects/Table_SNPs.Rds"))
  tab_length <- readRDS(file=file.path(path_to_output_folder,"00_R_objects/Table_consensus_length.Rds"))
}


tab_snps <- as.matrix(tab_snps)
loci <- t(tab_snps)

file.copy(from = config_file, to=file.path(output_assess,"0_used_config_file.txt"))

##########################################################
### Dataset optimization step 1: Reducing missing data ### 
##########################################################


nloci <- length(colnames(loci))
nsamples <- length(rownames(loci))

failed_loci <- which(colSums(is.na(loci))==nrow(loci))
failed_samples <- which(colSums(is.na(tab_snps))==nrow(tab_snps))

# per locus

seq_per_locus <- vector()
for(i in 1:nloci){
  seq_per_locus[i] <- length(which(!(is.na(loci[,i]))))
}
names(seq_per_locus) <- colnames(loci)
seq_per_locus_prop <- seq_per_locus/nsamples


# per sample

seq_per_sample <- vector()
for(i in 1:length(colnames(tab_snps))){
  seq_per_sample[i] <- length(which(!(is.na(tab_snps[,i]))))
}
names(seq_per_sample) <- colnames(tab_snps)
seq_per_sample_prop <- seq_per_sample/nloci


# proportion of target sequence length

if(targets_file_format == "AA"){
  targets_length_all <- lengths(read.FASTA(fasta_file_with_targets, type = "AA"))*3
} else if(targets_file_format == "DNA"){
  targets_length_all <- lengths(read.FASTA(fasta_file_with_targets))
} else { 
  print("Warning! Target file type not set properly. Should be 'DNA' or 'AA'!")
}


gene_names <- unique(gsub(".*-","",gsub(" .*","",names(targets_length_all))))
max_target_length <- vector()

for(i in 1:length(gene_names)){
  max_target_length[i] <- max(targets_length_all[grep(paste("\\b",gene_names[i],"\\b",sep=""),names(targets_length_all))])
}

names(max_target_length) <- gene_names
comb_target_length <- sum(max_target_length)
comb_seq_length_samples <- colSums(tab_length, na.rm = T)
prop_target_length_per_sample <- comb_seq_length_samples/comb_target_length
mean_seq_length_loci <- rowMeans(tab_length, na.rm = T)
names(mean_seq_length_loci) <- gene_names
prop_target_length_per_locus <- mean_seq_length_loci/max_target_length



# application of thresholds

outsamples_missing_loci <- seq_per_sample_prop[which(seq_per_sample_prop < remove_samples_with_less_than_this_propotion_of_loci_recovered)]
outsamples_missing_target <- prop_target_length_per_sample[which(prop_target_length_per_sample < remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered)]
outsamples_missing <- unique(names(c(outsamples_missing_loci,outsamples_missing_target)))

outloci_missing_samples <- seq_per_locus_prop[which(seq_per_locus_prop < remove_loci_with_less_than_this_propotion_of_samples_recovered)]
outloci_missing_target <- prop_target_length_per_locus[which(prop_target_length_per_locus < remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered)]
outloci_missing <- unique(names(c(outloci_missing_samples,outloci_missing_target)))


# removing bad loci and samples from the table 

tab_snps_cl1 <- tab_snps

if(length(outsamples_missing) != 0){
  tab_snps_cl1 <- tab_snps_cl1[,-which(colnames(tab_snps) %in% outsamples_missing)]
}

if(length(outloci_missing) != 0){
  tab_snps_cl1 <- tab_snps_cl1[-which(rownames(tab_snps) %in% outloci_missing),]
}

loci_cl1 <- t(tab_snps_cl1)



# output

# graphics

for(i in 1:2){
  
  if(i==1){
    pdf(file=file.path(output_assess,"1_Data_recovered_overview.pdf"), width = 11, height=7)
  } else {
    png(file=file.path(output_assess,"1_Data_recovered_overview.png"), width = 1400, height=1000)
    par(cex.axis=2, cex.lab=2, cex.main=2)
  }
  par(mfrow=c(2,3))
  
  boxplot(seq_per_sample_prop, main=paste("Samples: prop. of",nloci,"loci recovered"), xlab=paste("mean:",round(mean(seq_per_sample_prop, na.rm = TRUE),2), " | median:", round(median(seq_per_sample_prop, na.rm = TRUE),2)," | threshold:",remove_samples_with_less_than_this_propotion_of_loci_recovered," (",length(outsamples_missing_loci)," out)",sep=""))
  abline(h=remove_samples_with_less_than_this_propotion_of_loci_recovered, lty=2, col="red")
  
  boxplot(prop_target_length_per_sample, main=paste("Samples: prop. of target sequence length recovered"), xlab=paste("mean:",round(mean(prop_target_length_per_sample, na.rm = TRUE),2)," | median:", round(median(prop_target_length_per_sample, na.rm = TRUE),2)," | threshold:",remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered," (",length(outsamples_missing_target)," out)",sep=""))
  abline(h=remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  plot(prop_target_length_per_sample,seq_per_sample_prop, main = "Prop. of loci vs\n prop. of target length", xlab = "Prop. of target length", ylab= "Prop. of loci" )
  abline(h=remove_samples_with_less_than_this_propotion_of_loci_recovered, lty=2, col="red")
  abline(v=remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  
  boxplot(seq_per_locus_prop, main=paste("Loci: prop. of",nsamples,"samples recovered"), xlab=paste("mean:",round(mean(seq_per_locus_prop, na.rm = TRUE),2), " | median:", round(median(seq_per_locus_prop, na.rm = TRUE),2)," | threshold:",remove_loci_with_less_than_this_propotion_of_samples_recovered," (",length(outloci_missing_samples)," out)",sep=""))
  abline(h=remove_loci_with_less_than_this_propotion_of_samples_recovered, lty=2, col="red")
  
  boxplot(prop_target_length_per_locus, main=paste("Loci: prop. of target sequence length recovered"), xlab=paste("mean:",round(mean(prop_target_length_per_locus, na.rm = TRUE),2)," | median:", round(median(prop_target_length_per_locus, na.rm = TRUE),2)," | threshold:",remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered," (",length(outloci_missing_target)," out)",sep=""))
  abline(h=remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  plot(prop_target_length_per_locus,seq_per_locus_prop, main = "Prop. of samples vs\n prop. of target length", xlab = "Prop. of target length", ylab= "Prop. of samples" )
  abline(h=remove_loci_with_less_than_this_propotion_of_samples_recovered, lty=2, col="red")
  abline(v=remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  dev.off()
}


# tables

tab_seq_per_sample <- cbind(seq_per_sample,round(seq_per_sample_prop,3), round(prop_target_length_per_sample,3))
colnames(tab_seq_per_sample) <- c("No. loci", "Prop. of loci", "Prop. of target length")
write.csv(tab_seq_per_sample, file.path(output_assess, "1_Data_recovered_per_sample.csv"))

tab_seq_per_locus <- cbind(seq_per_locus,round(seq_per_locus_prop,3), round(prop_target_length_per_locus,3))
colnames(tab_seq_per_locus) <- c("No. samples", "Prop.of samples", "Prop. of target length")
write.csv(tab_seq_per_locus, file.path(output_assess, "1_Data_recovered_per_locus.csv"))


# summary text file
summary_file=file.path(output_assess,"1_Summary_missing_data.txt")
cat(file=summary_file, append = FALSE, "Dataset optimisation: Samples and loci removed to reduce missing data\n")

cat(file=summary_file, append = T, "\n", length(failed_samples)," samples failed completely:\n", paste(names(failed_samples)),"\n", sep="")
cat(file=summary_file, append = T, "\n", length(outsamples_missing_loci)," samples are below the threshold (",remove_samples_with_less_than_this_propotion_of_loci_recovered,") for proportion of recovered loci:\n", paste(names(outsamples_missing_loci),"\t",round(outsamples_missing_loci,3),"\n"), sep="")
cat(file=summary_file, append = T, "\n", length(outsamples_missing_target)," samples are below the threshold (",remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered,") for recovered target sequence length\n", paste(names(outsamples_missing_target),"\t",round(outsamples_missing_target,3),"\n"), sep="")
cat(file=summary_file, append = T, "\nIn total ", length(outsamples_missing), " samples were removed:\n", paste(outsamples_missing,"\n"), sep="")

cat(file=summary_file, append = T, "\n", length(failed_loci)," loci failed completely:\n", paste(names(failed_loci)),"\n", sep="")
cat(file=summary_file, append = T, "\n", length(outloci_missing_samples)," loci are below the threshold (", remove_loci_with_less_than_this_propotion_of_samples_recovered,") for proportion of recovered samples:\n", paste(names(outloci_missing_samples),"\t",round(outloci_missing_samples,3),"\n"), sep="")
cat(file=summary_file, append = T, "\n", length(outloci_missing_target)," loci are below the threshold (",remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered,") for proportion of recovered target sequence length:\n", paste(names(outloci_missing_target),"\t",round(outloci_missing_target,3),"\n"), sep="")
cat(file=summary_file, append = T, "\nIn total ", length(outloci_missing), " loci were removed:\n", paste(outloci_missing,"\n"), sep="")



############################################################################################
### Dataset optimization step 2, removing paralogs for a) all samples and b) each sample ### 
############################################################################################


### 2a) Paralogs across multiple samples (removing loci with unusually high proportions of SNPs across all samples)
###################################################################################################################

loci_cl1_colmeans <- colMeans(as.matrix(loci_cl1), na.rm = T)
nloci_cl1 <- length(colnames(loci_cl1))
nsamples_cl1 <- length(colnames(tab_snps_cl1))

loci_cl1_colmeans_mean <- round(mean(loci_cl1_colmeans),4)
loci_cl1_colmeans_median <- round(median(loci_cl1_colmeans),4)


# applying chosen threshold
if (length(remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs) != 1 || remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs == "none"){
  outloci_para_all <- vector()
  threshold_value <- 1
} else if (remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs == "outliers"){
  threshold_value <- 1.5*IQR(loci_cl1_colmeans, na.rm = TRUE )+quantile(loci_cl1_colmeans, na.rm = TRUE )[4]
  outloci_para_all_values <- loci_cl1_colmeans[which(loci_cl1_colmeans > threshold_value)]
  outloci_para_all <- names(outloci_para_all_values)
} else if (remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs == "file"){
  if(file.exists(file_with_putative_paralogs_to_remove_for_all_samples) == FALSE){
    print("File with list of paralogs to remove for all samples does not exist.")
  } else {
    threshold_value <- 1
    outloci_para_all <- readLines(file_with_putative_paralogs_to_remove_for_all_samples)
    outloci_para_all_values <-loci_cl1_colmeans[which(names(loci_cl1_colmeans) %in% outloci_para_all)]
  }
} else {
  threshold_value <- remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs
  outloci_para_all_values <- loci_cl1_colmeans[which(loci_cl1_colmeans > threshold_value)]
  outloci_para_all <- names(outloci_para_all_values)
}


# color outliers red
colour_outparaall <- rep("black",nloci_cl1)
colour_outparaall[which(colnames(loci_cl1[,order(loci_cl1_colmeans)]) %in% outloci_para_all)] <- "red"
loci_cl1_order_means <- loci_cl1[,order(loci_cl1_colmeans)]


# generate bar graph
for(i in 1:2){
  
  if(i==1){
    pdf(file=file.path(output_assess,"2a_Paralogs_for_all_samples.pdf"), width = 11, height=7)
  } else {
    png(file=file.path(output_assess,"2a_Paralogs_for_all_samples.png"), width = 1400, height=1000)
    par(cex.axis=2, cex.lab=2, cex.main=2)
  }
  
  layout(matrix(c(1,2),2,2, byrow=TRUE), widths=c(5,1))
  
  barplot(sort(loci_cl1_colmeans), col=colour_outparaall, border = NA, las=2,
          main=paste("Mean % SNPs across samples (n=",nsamples_cl1,") for each locus (n=", nloci_cl1,")", sep=""))
  if(length(threshold_value)>0 && remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs != "file"){abline(h=threshold_value, col="red", lty=2)}
  boxplot(loci_cl1_colmeans, las=2)
  if(length(threshold_value)>0 && remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs != "file"){abline(h=threshold_value, col="red", lty=2)}
  
  dev.off()
}


# removing marked loci from table
if(length(outloci_para_all)==0) {tab_snps_cl2a <- tab_snps_cl1
} else { tab_snps_cl2a <- tab_snps_cl1[-which(rownames(tab_snps_cl1) %in%  outloci_para_all),]}




### 2b) Paralogs for each sample (removing outlier loci for each sample)
##########################################################################


tab_snps_cl2b <- tab_snps_cl2a

if(!exists("remove_outlier_loci_for_each_sample")) {remove_outlier_loci_for_each_sample <- "no"}

# generate tables without zeros to count only loci with SNPs
tab_snps_cl2a_nozero <- tab_snps_cl2a
tab_snps_cl2a_nozero[which(tab_snps_cl2a_nozero==0)] <- NA
tab_snps_cl2b_nozero <- tab_snps_cl2a_nozero

if(remove_outlier_loci_for_each_sample == "yes" ){

  outloci_para_each <- list()
  outloci_para_each <- sapply(colnames(tab_snps_cl2a),function(x) NULL)
  threshold_para_each <- outloci_para_each
  for(i in 1:length(colnames(tab_snps_cl2a))){
    threshold_i <- 1.5*IQR(tab_snps_cl2a_nozero[,i], na.rm = TRUE ) + quantile(tab_snps_cl2a_nozero[,i] , na.rm = TRUE)[[4]]
    outlier_loci_i <- tab_snps_cl2a_nozero[which(tab_snps_cl2a_nozero[,i] > threshold_i),i]
    outloci_para_each[i] <- list(outlier_loci_i)
    threshold_para_each[i] <- threshold_i
    
    tab_snps_cl2b[which(rownames(tab_snps_cl2a) %in% outlier_loci_i),i] <- NA
    tab_snps_cl2b_nozero[which(rownames(tab_snps_cl2b_nozero) %in% names(outlier_loci_i)),i] <- NA
    
    outliers_color <- "red"
  }
  
} else {
  outloci_para_each <- list()
  outloci_para_each <- sapply(colnames(tab_snps_cl2a),function(x) NULL)
  tab_snps_cl2b <- tab_snps_cl2a
  outliers_color <- "black"
}

tab_length_cl2b <- tab_length[which(rownames(tab_length) %in% rownames(tab_snps_cl2b)),which(colnames(tab_length) %in% colnames(tab_snps_cl2b))]


### output 

# generate graphic
for (i in 1:2){
  if(i==1){
    pdf(file=file.path(output_assess,"2b_Paralogs_for_each_sample.pdf"), width = 10, h=14)
  } else {
    png(file=file.path(output_assess,"2b_Paralogs_for_each_sample.png"), width = 1000, h=1400)
  }
  par(mfrow=c(1,1))
  boxplot(as.data.frame(tab_snps_cl2a_nozero[,order(colMeans(as.matrix(tab_snps_cl2a_nozero), na.rm = T))]), 
          horizontal=T, las=1, yaxt='n',
          main="Proportions of SNPs for all loci per sample\n(only loci with any SNPs)",
          ylab="Samples",
          xlab="Proportion of SNPs",
          col="grey", pars=list(outcol=outliers_color),
          outpch=20
          )
  dev.off()  
}



# write summary text file 
cl2b_file <- file.path(output_assess,"2_Summary_Paralogs.txt")
cat(file=cl2b_file,"Removal of putative paralog loci.")
cat(file=cl2b_file,"Paralogs removed for all samples:\n", append = T)
cat(file=cl2b_file, paste("Variable 'remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs' set to: ", remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs,"\n", sep=""), append=T)
if(remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs=="file"){
  cat(file=cl2b_file, paste("Loci listed in this file were removed: '", file_with_putative_paralogs_to_remove_for_all_samples,"'\n"), append=T)
} else if(remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs=="none"){
  cat(file=cl2b_file,"None!\n", append = T)
} else {
  cat(file=cl2b_file, paste("Resulting threshold value (mean proportion of SNPs):", round(threshold_value,5),"\n"), append=T)
}
if(length(outloci_para_all)>0){
  cat(file=cl2b_file, paste(length(outloci_para_all)," loci were removed:\n",sep=""), append=T)
  cat(file=cl2b_file, "locus\tmean_prop_SNPs\n", append=T)
  cat(file=cl2b_file, paste(paste(outloci_para_all, round(outloci_para_all_values,4), sep="\t"), collapse="\n"), append=T)
  cat(file=cl2b_file,  "\n\n", append = T)
}

cat(file= file.path(output_assess,"2a_List_of_paralogs_removed_for_all_samples.txt"), paste(c(outloci_para_all,"\n"),collapse = "\n"))



if(remove_outlier_loci_for_each_sample=="yes"){
  cat(file=cl2b_file, "Paralogs removed for each sample:\n", append=T)
  cat(file=cl2b_file, "Sample\tthreshold\t#removed\tnames\n", append=T)
  for(i in 1:length(names(outloci_para_each))){
    cat(file=cl2b_file, names(outloci_para_each)[i],"\t", append=T)  
    cat(file=cl2b_file, round(threshold_para_each[[i]],5), length(outloci_para_each[[i]]), paste(names(outloci_para_each[[i]]),collapse=", "), sep="\t", append=T)  
    cat(file=cl2b_file, "\n", append=T)  
  }
} else{
  cat(file=cl2b_file, "The step for removing paralogs for each samples was skipped.\n", append=T) 
}


# tables

write.csv(tab_snps_cl2b, file = file.path(output_assess,"0_Table_SNPs.csv"))
write.csv(tab_length_cl2b, file = file.path(output_assess,"0_Table_consensus_length.csv"))

# txt file with included samples
write(rownames(loci)[which(!(rownames(loci) %in% outsamples_missing))], file=file.path(output_assess,"0_namelist_included_samples.txt"))

#write(paste(rownames(tab_snps_cl2a),colMeans(t(tab_snps_cl2a), na.rm = T)), file=file.path(output_assess,"Mean_SNPs_loci.txt"))


### save Data as R objects

saveRDS(tab_snps_cl2b,file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))
saveRDS(tab_length_cl2b,file=file.path(output_Robjects,"Table_consensus_length_cleaned.Rds"))

saveRDS(outloci_missing,file=file.path(output_Robjects,"outloci_missing.Rds"))
saveRDS(outsamples_missing,file=file.path(output_Robjects,"outsamples_missing.Rds"))
saveRDS(outloci_para_all,file=file.path(output_Robjects,"outloci_para_all.Rds"))
saveRDS(outloci_para_each,file=file.path(output_Robjects,"outloci_para_each.Rds"))



#############################################################################################################
### Generating summary table and graphs for assessment of Locus heterozygosity and allele divergence of samples ###
#############################################################################################################


tab_length <- as.matrix(tab_length)

tab_length_cl2b <- tab_length[which(rownames(tab_length) %in% rownames(tab_snps_cl2b)),which(colnames(tab_length) %in% colnames(tab_snps_cl2b))]

targets_length_cl2b <- sum(max_target_length[which(gsub(".*-","",names(max_target_length)) %in% rownames(tab_snps_cl2b))])

########## generating summary table 

nloci_cl2 <- length(tab_snps_cl2b[,1])
tab_het_ad <- data.frame("sample"=colnames(tab_snps_cl2b))

for(i in 1:length(colnames(tab_snps_cl2b))){
  tab_het_ad$bp[i] <- sum(tab_length_cl2b[,i], na.rm = T)
  tab_het_ad$bpoftarget[i] <- round(sum(tab_length_cl2b[,i], na.rm = T)/targets_length_cl2b,3)*100
  tab_het_ad$paralogs_all[i] <- length(outloci_para_all)
  tab_het_ad$paralogs_each[i] <- length(outloci_para_each[[i]])
  tab_het_ad$nloci[i] <- nloci_cl1-length(outloci_para_all)-length(which(is.na(tab_snps_cl2b[,i])))-length(outloci_para_each[[i]])
  tab_het_ad$allele_divergence[i] <- 100*round(sum(tab_length_cl2b[,i] * tab_snps_cl2b[,i], na.rm = T) / sum(tab_length_cl2b[,i], na.rm = T),5)
  tab_het_ad$locus_heterozygosity[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]==0))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'loci with >0.5% SNPs'[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]<0.005))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'loci with >1% SNPs'[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]<0.01))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'loci with >2% SNPs'[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]<0.02))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
}

# output as csv file
write.csv(tab_het_ad, file = file.path(output_assess, "4_Summary_table.csv"))

#output as R-object
saveRDS(tab_het_ad, file = file.path(output_Robjects, "Summary_table.Rds"))



### Generating graphs


text_size_mod <- 1
nrows <- length(tab_het_ad[,1])
text_size <- (15+200/nrows)*text_size_mod


for(i in 1:2){
  
  if(i==1){
    pdf(file.path(output_assess,"3_LH_vs_AD.pdf"), h=10,w=10)
  } else {
    png(file.path(output_assess,"3_LH_vs_AD.png"), h=1000,w=1000)
  }
  
  plot(tab_het_ad$allele_divergence,tab_het_ad$locus_heterozygosity,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity [%]", main="Locus heterozygosity vs allele divergence", las=1)
  dev.off()
}


for(i in 1:2){
  
  if(i==1){
    pdf(file.path(output_assess,"3_varLH_vs_AD.pdf"), h=10,w=10)
  } else {
    png(file.path(output_assess,"3_varLH_vs_AD.png"), h=1000,w=1000)
  }
  
  par(mfrow=c(2,2))
  plot(tab_het_ad$allele_divergence,tab_het_ad$locus_heterozygosity,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (0% SNPs) [%]", main="Locus heterozygosity (0% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence,tab_het_ad$`loci with >0.5% SNPs`,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (>0.5% SNPs) [%]", main="Locus heterozygosity (>0.5% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence,tab_het_ad$`loci with >1% SNPs`,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (>1% SNPs) [%]", main="Locus heterozygosity (>1% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence,tab_het_ad$`loci with >2% SNPs`,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (>2% SNPs) [%]", main="Locus heterozygosity (>2% SNPs) vs allele divergence",las=1
  )
  par(mfrow=c(1,1))
  dev.off()
}




