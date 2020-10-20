# HybPhaser Main script

# This main script is used to read all configuration scripts and execute relevant subscripts. 
# Make sure you set all variables in the configuration scripts before running them. 


### Packages
library(ape)
library(seqinr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)


##################################
### Part 1: Assessment of SNPs ###
##################################

# Set all variables in the configuration scripts before running next line!
source("Config_1_SNPs_assessment.R")

### run the next line to execute script to count SNPs in consensus files (this will take a few minutes)
source("Rscript_1a_count_snps_in_consensus_seqs.R")

### run next line to excute the script for dataset optimization 
source("Rscript_1b_optimize_dataset.R")

# check output (".../HybPhser_folder/dataset_optimization/") and if required, adjust thresholds and rerun dataset optimization

### run next line to generate graphs and tables for assessment of heterozygosity and allele divergence
source("Rscript_1c_summary_table.R")

### run next line to run script to generate sequence lists 
source("Rscript_1d_generate_sequence_lists.R")

# Sequence lists are available in subfolder and ready for alignment and phylogeny reconstruction


##########################################
### Part 2: Clade Association Analysis ###
##########################################

# Set all variables in the configuration script before running next line!
source("Configure_2_Clade_association.R")


### run next line to prepare (and optionally run) clade association 
source("Rscript_2a_prepare_bbsplit_script.R")

# if the bash script was not run in R, execute the bash script in the command line before running the next line

## run next line after BBSplit has finished to collate results
source("Rscript_2b_collate_bbsplit_results.R")


# The tables with collated results are in the clade association folder and ready for review! 



#########################################
### Part 3: Phasing of sequence reads ###
#########################################

# Set all variables in the configuration script before running next line!
source("Configure_3_Phasing.R")


### run next line to prepare (and optionally run) phasing script
source("Rscript_3a_prepare_phasing_script.R")

# if the bash script was not run in R, execute the bash script in the command line to generate read files for the newly phased accessions

# run next line to collate phasing stats into a table
source("Rscript_3b_collate_phasing_stats.R")

# Tables with phasing stats are available for review 
# Phased sequence read files are avaiable and ready to run through HybPiper and HybPhaser


############################################################################
### Part 4: Combining Sequence Lists of phased and non-phased accessions ###
############################################################################


# set all variables in the configuration script before running next line!
source("Configure_4_Combining_phased_with_normal_sequence_lists.R")

### execute script to generate new combined sequence lists
source("Rscript_4a_combine_phased_with_normal_sequence_lists.R")

