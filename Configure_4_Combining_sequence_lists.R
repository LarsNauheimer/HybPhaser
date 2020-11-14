###################################################################
### Configuration for merging phased with normal sequence lists ###
####################################################################

# This script is to generate sequence lists that combine the phased sequences with the normal non-phased ones (but exclude the normal ones of the phased accessions)
# It is also possible to make subsets of samples or loci using lists of included or excluded samples/loci


# set HybPhaser base folder
path_for_HybPhaser_output = ""

# set HybPhaser base folder of phased accessions
path_for_HybPhaser_phased_output = ""

# set sequence list name (referring to the folder name in the sequences subfolder, e.g. "consensus_loci_clean" or "contig_loci_clean") 
sequence_type = ""

# set contig type ("normal" or "supercontig")
contig = ""    

# name the output sequence list folder 
name_of_sequence_list_output = ""

# set nameslist file of original accessions (e.g. "~/hybpiper/namelist.txt")
txt_file_with_list_of_accessions = ""

# set nameslist file of phased accessions (e.g. "~/hybpiper_phased/namelist_phased.txt")
txt_file_with_list_of_phased_accessions = ""

# select whether the phased accessions should be used instead of the non phased accessions of a sample. If "no" is chosen the phased and not phased accessions of the same sample are included. 
exchange_phased_with_not_phased_samples = "yes"   # "yes" or "no" 

# finally, sometimes it can be that there are phased loci included that were excluded in the non-phased list. In case they should be included, set the next variable to "yes"
include_phased_seqlists_when_non_phased_locus_absent = "no"  # "yes" or "no"


# Optional: adjust list of samples and loci 
#
# The default is to use the sample lists used in previous scripts, one for the normal and one for the phased reads. 
# They will be merged and all samples that appear in both lists will be used. 
# Samples can be removed by listing them in text file (file_with_samples_excluded).
#
# Or one can provide a list only with the included samples (file_with_samples_included).
# The other option is to provide a text file with included samples and thus overwriting the default. 
#
# The same can be done for loci.

# select samples to exclude  (one sample per row in a text file).  "" will exclude none.
file_with_samples_excluded = ""

# Or alternatively only use one file that lists all samples included (phased and non-phased).
file_with_samples_included = ""   


# Similarly one can choose to exclude loci or use a list with only the included loci. 

# select loci to exclude (text file with one locus per row). "" will exclude none.
file_with_loci_excluded = ""

# select loci to include (text file with one locus per row). "" will include all.
file_with_loci_included = ""

