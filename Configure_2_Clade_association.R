############################################
### Configuration for clade association ###
############################################

# set variables to determine thresholds for the dataset optimization and run the script from the main script

path_to_HybPhaser_results = ""
path_to_clade_association_folder = ""
csv_file_with_clade_reference_names = ""
path_to_read_files_cladeassociation = ""
read_type_cladeassociation = "single-end" 
ID_read_pair1 = ""
ID_read_pair2 = ""
file_with_samples_included <- ""
path_to_reference_sequences = ""
path_to_bbmap = ""
no_of_threads = 1
run_bash_script_in_R = "no"   


# Description fo variables
##########################

# path_to_HybPhaser_results
## path to HybPhaser result folder (r.g. .../HybPiper/_HybPhaser/)

# path_to_clade_association_folder
## set clade_association folder. e.g. as subfolder in the HybPhaser folder ("~/HybPhaser_supercontig/clade_association_1")

# csv_file_with_clade_reference_names
## set filename for the previously prepared text file with reference sample names (and optionally abbreviations for each reference)) 
## (e.g. "~/HybPhaser_supercontig/clade_association_1/clade_references.csv")

# path_to_read_files_cladeassociation
## set folder that contains read files (e.g. "~/HybPiper/_mapped_reads/")
## it is recommended to use only the reads mapped to the targets. They can be extracted using the bash script for extracting reads from the HybPiper bam file

# read_type_cladeassociation
## set whether reads are single or paired-end reads (single end, if you use the mapped-only reads

# ID_read_pair1 / ID_read_pair2
## if reads are paird-end, set unique part of filename including the ending (e.g. "_R1.fastq", "_R2.fastq")
## if reads are single-end, you can ignore these 

# file_with_samples_included
## set path to text file that contains a list of all samples included

# path_to_reference_sequences
## set folder that contains the consensus sequences of samples (e.g. ".../HybPiper_results/HybPhaser/sequences/consensus_samples_raw")

# path_to_bbmap
## set folder to bbmap binaries, if bbmap is not in your path variable (e.g. "~/applications/bbmap")

# no_of_threads
## set number of threads to be used

# run_bash_script_in_R
## select whether the script is run directly from R (it might take a while and might be better to be run in a separate terminal window)


