###################################################################
### Configuration for merging phased with normal sequence lists ###
####################################################################

# This script is to generate sequence lists that combine the phased sequences with the normal non-phased ones (but exclude the normal ones of the phased accessions)
# It is also possible to make subsets of samples or loci using lists of included or excluded samples/loci

path_for_HybPhaser_output = ""
path_for_HybPhaser_phased_output = ""
sequence_type = ""
contig = ""    
name_of_sequence_list_output = ""

# Sample list
# The default is to use the sample lists used in previous scripts, one for the normal and one for the phased reads. 
# They will be merged and all samples that appear in both lists as well as are available in the sequence lists will be used. 
# From that list samples can be removed by listing them in text file (samples to exclude).
txt_file_with_list_of_accessions = ""
txt_file_with_list_of_phased_accessions = ""
file_with_samples_excluded = ""

# The other option is to provide a text file with included samples and thus overwriting the default. 
file_with_samples_included = ""   

# Similarly one can choose to exclude loci or use a list with only the included loci. 
file_with_loci_excluded = ""
file_with_loci_included = ""

exchange_phased_with_not_phased_samples = "yes"   # "yes" or "no" 
include_phased_seqlists_when_non_phased_locus_absent = "no"  # "yes" or "no"


# Description fo variables
##########################

# # path_for_HybPhaser_output
## set HybPhaser base folder

# path_for_HybPhaser_phased_output
## set HybPhaser base folder of phased accessions

# sequence_type
## set sequence list name (referring to the folder name in the sequences subfolder, e.g. "consensus_loci_clean" or "contig_loci_clean") 

#  contig
## set contig type ("normal" or "supercontig")

# name_of_sequence_list_output
## name the output sequence list folder 

# Sample list
# The default is to use the sample lists used in previous scripts, one for the normal and one for the phased reads. 
# They will be merged and all samples that appear in both lists as well as are available in the sequence lists will be used. 
# From that list samples can be removed by listing them in text file (samples to exclude).
# The other option is to provide a text file with included samples and thus overwriting the default. 

# txt_file_with_list_of_accessions
## set nameslist file of original accessions (e.g. "~/hybpiper/namelist.txt")

# txt_file_with_list_of_phased_accessions
## set nameslist file of phased accessions (e.g. "~/hybpiper_phased/namelist_phased.txt")

# file_with_samples_excluded
## select samples to exclude  (one sample per row in a text file).  "" will exclude none.

# file_with_samples_included
## Or alternatively only use one file that lists all samples included (phased and non-phased).

# file_with_loci_excluded
## select loci to exclude (text file with one locus per row). "" will exclude none.

# file_with_loci_included
## select loci to include (text file with one locus per row). "" will include all.

# exchange_phased_with_not_phased_samples
## select whether the phased accessions should be used instead of the non phased accessions of a sample. If "no" is chosen the phased and not phased accessions of the same sample are included. 

# include_phased_seqlists_when_non_phased_locus_absent
## finally, sometimes it can be that there are phased loci included that were excluded in the non-phased list. In case they should be included, set the next variable to "yes"

