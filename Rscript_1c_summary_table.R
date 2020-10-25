#############################################################################################################
### Generating summary table and graphs for assessment of heterozygosity and allele divergence of samples ###
#############################################################################################################

#input
tab_snps_cl2b <- readRDS(file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))
tab_length <- readRDS(file=file.path(output_Robjects,"Table_consensus_length.Rds"))
tab_length_cl2b <- readRDS(file=file.path(output_Robjects,"Table_consensus_length_cleaned.Rds"))

outloci_para_all <- readRDS(file=file.path(output_Robjects,"outloci_para_all.Rds"))
outloci_para_each <- readRDS(file=file.path(output_Robjects,"outloci_para_each.Rds"))

dir.create(output_assess, showWarnings = F)

rownames(tab_length) <-  tab_length[,1]
tab_length <- as.matrix(tab_length[,-c(1)])

tab_length_cl2b <- tab_length[which(rownames(tab_length) %in% rownames(tab_snps_cl2b)),which(colnames(tab_length) %in% colnames(tab_snps_cl2b))]

targets_length_all <- lengths(read.FASTA(fasta_file_with_targets))
gene_names <- unique(gsub(".*-","",names(targets_length_all)))
max_target_length <- vector()

for(i in 1:length(gene_names)){
  max_target_length[i] <- max(targets_length_all[grep(gene_names[i],names(targets_length_all))])
}

names(max_target_length) <- gene_names

targets_length <- sum(max_target_length)
targets_length_cl2b <- sum(max_target_length[which(gsub(".*-","",names(max_target_length)) %in% rownames(tab_snps_cl2b))])

########## generating summary table 

nloci_cl2 <- length(tab_snps_cl2b[,1])
tab_het_ad <- data.frame("samples"=colnames(tab_snps_cl2b))

for(i in 1:length(colnames(tab_snps_cl2b))){
  tab_het_ad$bp[i] <- sum(tab_length[,grep(colnames(tab_snps_cl2b)[i],colnames(tab_length))], na.rm = T)
  tab_het_ad$bpoftarget[i] <- round(sum(tab_length[,grep(colnames(tab_snps_cl2b)[i],colnames(tab_length))], na.rm = T)/targets_length,3)*100
  tab_het_ad$bp_clean[i] <- sum(tab_length_cl2b[,i], na.rm = T)
  tab_het_ad$bpoftarget_clean[i] <- round(sum(tab_length_cl2b[,i], na.rm = T)/targets_length_cl2b,3)*100
  tab_het_ad$paralogs_all[i] <- length(outloci_para_all)
  tab_het_ad$paralogs_each[i] <- length(outloci_para_each[[i]])
  tab_het_ad$nloci[i] <- nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))
  tab_het_ad$'>0% SNPs'[i] <- round(1 - length(which(tab_snps_cl2b[,i]==0))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'>0.1% SNPs'[i] <- round(1 - length(which(tab_snps_cl2b[,i]<0.001))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'>0.5% SNPs'[i] <- round(1 - length(which(tab_snps_cl2b[,i]<0.005))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'>1% SNPs'[i] <- round(1 - length(which(tab_snps_cl2b[,i]<0.01))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'>2% SNPs'[i] <- round(1 - length(which(tab_snps_cl2b[,i]<0.02))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$allele_divergence[i] <- round(sum(tab_length_cl2b[,i] * tab_snps_cl2b[,i], na.rm = T) / sum(tab_length_cl2b[,i], na.rm = T),4)
}


# output as csv file
write.csv(tab_het_ad, file = file.path(output_assess, "Table_Summary.csv"))

#output as R-object
saveRDS(tab_het_ad, file = file.path(output_Robjects, "Table_Summary.Rds"))


### Generating graphs


text_size_mod <- 1
nrows <- length(tab_het_ad[,1])
text_size <- (15+200/nrows)*text_size_mod


pdf(file.path(output_assess,"Scatterplot_heterozygosity_vs_allele_divergence.pdf"), h=10,w=10)
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>0% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity [%]", main="Heterozygosity vs allele divergence", las=1)
dev.off()

png(file.path(output_assess,"Scatterplot_heterozygosity_vs_allele_divergence.png"), h=1000,w=1000)
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>0% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity [%]", main="Heterozygosity vs allele divergence",las=1)
dev.off()


pdf(file.path(output_assess,"Scatterplot_heterozygosity_div_levels_vs_allele_divergence.pdf"), h=10,w=10)
  par(mfrow=c(2,2))
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>0% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (0% SNPs) [%]", main="Heterozygosity (0% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>0.5% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (>0.5% SNPs) [%]", main="Heterozygosity (.0.5% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>1% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (>1% SNPs) [%]", main="Heterozygosity (>1% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>2% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (>2% SNPs) [%]", main="Heterozygosity (>2% SNPs) vs allele divergence",las=1
  )
  par(mfrow=c(1,1))
dev.off()           


png(file.path(output_assess,"Scatterplot_heterozygosity_div_levels_vs_allele_divergence.png"), h=1000,w=1000)
  par(mfrow=c(2,2))
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>0% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (0% SNPs) [%]", main="Heterozygosity (0% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>0.5% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (>0.5% SNPs) [%]", main="Heterozygosity (.0.5% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>1% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (>1% SNPs) [%]", main="Heterozygosity (>1% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence*100,tab_het_ad$`>2% SNPs`*100,
       xlab="Allele divergence [%]", ylab="Heterozygosity (>2% SNPs) [%]", main="Heterozygosity (>2% SNPs) vs allele divergence",las=1
  )
  par(mfrow=c(1,1))
dev.off()

