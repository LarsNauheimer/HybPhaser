# HybPhaser Main script

# This main script is used to read all configuration scripts and execute relevant subscripts. 
# Make sure you set all variables in the configuration scripts before running them. 

### Packages ###
library(ape)
library(seqinr)
library(stringr)
################

path_to_config_scripts <- "."
path_to_Rscripts <- "."

##################################
### Part 1: Assessment of SNPs ###
##################################

### Configure Part 1
# Set all variables in the configuration scripts before running next line!
source(file.path(path_to_config_scripts,"Configure_1_SNPs_assessment.R"))


### Execute Part 1

# 1a) execute script to count SNPs in consensus files (this will take a few minutes)
source(file.path(path_to_Rscripts,"Rscript_1a_count_snps_in_consensus_seqs.R"))

# 1b) excute the script for dataset optimization 
source(file.path(path_to_Rscripts,"Rscript_1b_optimize_dataset.R"))

# check output (".../HybPhser_folder/dataset_optimization/") and if required, adjust thresholds and rerun dataset optimization

# 1c) generate graphs and tables for assessment of heterozygosity and allele divergence
source(file.path(path_to_Rscripts,"Rscript_1c_summary_table.R"))

# 1d) generate sequence lists 
source(file.path(path_to_Rscripts,"Rscript_1d_generate_sequence_lists.R"))

# Sequence lists are available in subfolder and ready for alignment and phylogeny reconstruction


##########################################
### Part 2: Clade Association Analysis ###
##########################################

### Configure Part 2

# Set all variables in the configuration script before running next line!
source(file.path(path_to_config_scripts,"Configure_2_Clade_association.R"))


### Execute Part 2

# prepare (and optionally run) clade association 
source(file.path(path_to_Rscripts,"Rscript_2a_prepare_bbsplit_script.R"))

# if the bash script was not run in R, execute the bash script in the command line before running the next line

## run next line after BBSplit has finished to collate results
source(file.path(path_to_Rscripts,"Rscript_2b_collate_bbsplit_results.R"))


# The tables with collated results are in the clade association folder and ready for review! 



#########################################
### Part 3: Phasing of sequence reads ###
#########################################

### Configure Part 3

# Set all variables in the configuration script before running next line!
source(file.path(path_to_config_scripts,"Configure_3_Phasing.R"))


### Execute Part 3

### prepare (and optionally run) phasing script
source(file.path(path_to_Rscripts,"Rscript_3a_prepare_phasing_script.R"))

# if the bash script was not run in R, execute the bash script in the command line to generate read files for the newly phased accessions

# collate phasing stats into a table
source(file.path(path_to_Rscripts,"Rscript_3b_collate_phasing_stats.R"))

# Tables with phasing stats are available for review 
# Phased sequence read files are avaiable and ready to run through HybPiper and HybPhaser


############################################################################
### Part 4: Combining Sequence Lists of phased and non-phased accessions ###
############################################################################

### Configure Part 4

# set all variables in the configuration script before running next line!
source(file.path(path_to_config_scripts,"Configure_4_Combining_phased_with_normal_sequence_lists.R"))


### Execute Part 4

### generate new combined sequence lists
source(file.path(path_to_Rscripts,"Rscript_4a_combine_phased_with_normal_sequence_lists.R"))

