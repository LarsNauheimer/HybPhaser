############################################
### Configuration for clade association ###
############################################

# set variables to determine thresholds for the dataset optimization and run the script from the main script

# path to HybPhaser result folder (r.g. .../HybPiper/_HybPhaser/)
path_to_HybPhaser_results = ""

# set clade_association folder. e.g. as subfolder in the HybPhaser folder ("~/HybPhaser_supercontig/clade_association_1")
path_to_clade_association_folder = ""

# set filename for the previously prepared text file with reference sample names (and optionally abbreviations for each reference)) 
# (e.g. "~/HybPhaser_supercontig/clade_association_1/clade_references.csv")
csv_file_with_clade_reference_names = ""

# set folder that contains read files (e.g. "~/HybPiper/_mapped_reads/")
# it is recommended to use only the reads mapped to the targets. They can be extracted using the bash script for extracting reads from the HybPiper bam file
path_to_read_files_cladeassociation = ""

# set whether reads are single or paired-end reads (single end, if you use the mapped-only reads
read_type_cladeassociation = "single-end" # "single-end" or "paired-end"

# if reads are paird-end, set unique part of filename including the ending (e.g. "_R1.fastq", "_R2.fastq")
# if reads are single-end, you can ignore this 
ID_read_pair1 = "_1.fastq.gz"
ID_read_pair2 = "_2.fastq.gz"

file_with_samples_included <- ""

# set folder that contains the consensus sequences of samples (e.g. ".../HybPiper_results/HybPhaser/sequences/consensus_samples_raw")
path_to_reference_sequences = ""

# set folder to bbmap binaries, if bbmap is not in your path variable (e.g. "~/applications/bbmap")
path_to_bbmap = ""

# set number of threads to be used
no_of_threads = 1

# select whether the script is run directly from R (it might take a while and might be better to be run in a separate terminal window)
run_bash_script_in_R = "no"   # "yes" or "no"

