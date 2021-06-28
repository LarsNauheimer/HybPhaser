# HybPhaser Main script to run all other R scripts in R-Studio

# Make sure you set all variables in the config file before running the scripts. 



# set path to config file
config_file="config.txt"


##################################
### Part 1: Assessment of SNPs ###
##################################


# 1a) execute script to count SNPs in consensus files (this will take a few minutes)
source("1a_count_snps.R")

# 1b) excute the script to generate graphs for the data assessment 
source("1b_assess_dataset.R")

# check output ("output_folder/assessment/") and if required, adjust thresholds and rerun script 1b (maybe under different subset name)

# 1c) generate sequence lists 
source("1c_generate_sequence_lists.R")

# Sequence lists are available in subfolder and ready for alignment and phylogeny reconstruction


##########################################
### Part 2: Clade Association Analysis ###
##########################################

# prepare (and optionally run) clade association 
source("2a_prepare_bbsplit_script.R")

# if the bash script was not run in R, execute the bash script in the command line before running the next line

## run next line after BBSplit has finished to collate results
source("2b_collate_bbsplit_results.R")


# The tables with collated results are in the clade association folder and ready for review! 



#########################################
### Part 3: Phasing of sequence reads ###
#########################################

### prepare (and optionally run) phasing script
source("3a_prepare_phasing_script.R")

# if the bash script was not run in R, execute the bash script in the command line to generate read files for the newly phased accessions

# collate phasing stats into a table
source("3b_collate_phasing_stats.R")

# Tables with phasing stats are available for review 
# Phased sequence read files are avaiable and ready to run through HybPiper and HybPhaser


############################################################################
### Part 4: Combining Sequence Lists of phased and non-phased accessions ###
############################################################################

### generate new combined sequence lists
source("4_merge_sequence_lists.R")

