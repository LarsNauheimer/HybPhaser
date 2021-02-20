#######################################
### HybPhaser configuration for phasing
#######################################
# set variables and direct the scripts to the right folders and files before running the single scripts (from inside this file)

path_to_HybPhaser_results = ""
path_to_read_files_phasing = ""
read_type_4phasing = "single-end"   
ID_read_pair1 = ""
ID_read_pair2 = ""
reference_sequence_folder = "consensus_samples_raw"
path_to_bbmap_executables = ""
path_to_phasing_folder = ""
csv_file_with_phasing_prep_info = ""
no_of_threads_phasing = 1 
run_bash_script_in_R = "yes"      
folder_for_phased_reads = ""
folder_for_phasing_stats = ""


# Description fo variables
##########################

# path_to_HybPhaser_results
## set HybPhaser basefolder

# path_to_read_files_phasing
## set path to read files that should be phased

# read_type_4phasing
## set whether reads are paired-end or single-end
 
# ID_read_pair1 / ID_read_pair2
## if reads are paird-end, set unique part of filename including the ending (e.g. "_R1.fastq", "_R2.fastq")
 #reference_sequence_folder
## set folder for reference sequences (e.g. "consensus_samples_raw"). It must be one of the generated sequence lists for samples.
 
# path_to_bbmap_executables
## set path to bbmap executables (if not in path)
 
# path_to_phasing_folder
## set path to folder with phaseing results (e.g. "~/HybPhaser/phasing)

# csv_file_with_phasing_prep_info
## set path to preparation file (e.g. "~/HybPhaser/phasing/phasing_prep.csv")

# no_of_threads_phasing
## set number of maximum used threads 

# run_bash_script_in_R
## set whether to run the script to run the bbsplit analyses in R instead of in the terminal

# folder_for_phased_reads
## output folder for phased reads and phasing stats (e.g. "~/HybPhaser/phasing/phased_reads")

# folder_for_phasing_stats
## output folder for phasing stats (e.g. "~/HybPhaser/phasing/phasing_stats")


