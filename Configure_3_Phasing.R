#######################################
### HybPhaser configuration for phasing
#######################################
# set variables and direct the scripts to the right folders and files before running the single scripts (from inside this file)
 
# set HybPhaser basefolder
path_to_HybPhaser_results = ""

# set path to read files that should be phased
path_to_read_files_phasing = ""

# set whether reads are paired-end or single-end
read_type_4phasing = "single-end"       # "single-end" or "paired-end"

# if reads are paird-end, set unique part of filename including the ending (e.g. "_R1.fastq", "_R2.fastq")
ID_read_pair1 = "_R1.fastq"
ID_read_pair2 = "_R2.fastq"

# set folder for reference sequences (e.g. "consensus_samples_raw"). It must be one of the generated sequence lists for samples.
reference_sequence_folder = "consensus_samples_raw"

# set path to bbmap executables (if not in path)
path_to_bbmap_executables = ""

# set path to folder with phaseing results (e.g. "~/HybPhaser/phasing)
path_to_phasing_folder = ""


# prepare csv file with phasable samples and the references it should be phased to:
# comma separated file with the following header (columns):
# samples,ref1,abb1,ref2,abb2,ref3,abb3,...
# samples is the sample name that should be phased
# ref1 is the sample name of the 1st reference
# abb1 is the abbreviation for the 1st reference
# ...
# 
# set path to preparation file (e.g. "~/HybPhaser/phasing/phasing_prep.csv")
csv_file_with_phasing_prep_info = ""

# set number of maximum used threads 
no_of_threads_phasing=4

#set whether to run the script to run the bbsplit analyses in R instead of in the terminal
run_bash_script_in_R="yes"       # "yes" or "no"

# output folder for phased reads and phasing stats (e.g. "~/HybPhaser/phasing/phased_reads")
folder_for_phased_reads = ""

# output folder for phasing stats (e.g. "~/HybPhaser/phasing/phasing_stats")
folder_for_phasing_stats = ""

