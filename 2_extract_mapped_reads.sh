#!/bin/bash
############
set +e
set -u
set -o pipefail



BASEDIR="./"
OUTDIR="./mapped_reads/"
DEDUPSEQ="FALSE"
NAMELIST="not-a-file"

while getopts 'b:o:sn:' OPTION; do
  case "$OPTION" in
   
    b)
      BASEDIR=$OPTARG
     ;;         

    o)
      OUTDIR=$OPTARG
     ;;         

    s)
      DEDUPSEQ="TRUE"
     ;;         

    n)
     NAMELIST=$OPTARG
     ;;         
 

    ?)
      echo "script usage: $(basename $0) [-n file with samples names (nameslist.txt)] [-b HybPiper base folder] [-o output folder]"
      echo " 
      
      Usage: 2_extract_mapped_reads.sh [options]
      
        Options:
    
        -b <PATH> Base folder of either HybPhaser (contains '01_data/'), HybPiper (contains sample directories), or HybPiper-RBGV (contains '04_processed_gene_directories/'). Default is './'
            
        -o <PATH> Outpout folder to write read files into. Default is './mapped_reads/'
        
        -n <PATH> Namelist. optional, if not set, all folders in samples directory are used as samples. 
                
        -s Select if duplicated sequences should be removed. "
     
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"


mkdir -p $OUTDIR

if [[ -d $BASEDIR/01_data ]]; then 

	SAMPLESDIR=$BASEDIR/01_data/
	HYB="phaser"
	echo $HYB

elif [[ -d $BASEDIR/04_processed_gene_directories/ ]]; then 
	
	SAMPLESDIR=$BASEDIR/04_processed_gene_directories/
	HYB="piper"
	
else 
	
	SAMPLESDIR=$BASEDIR
	HYB="piper"
	
fi



if [[ -f $NAMELIST ]]; then
	
	SAMPLES=$(<$NAMELIST)

else
	
	SAMPLES=$(ls -d -1 $SAMPLESDIR/* | rev | cut -f 1 -d "/" | rev)

fi


for SAMPLE in $SAMPLES
do
	echo Generating on-target reads files for $SAMPLE

	# concatenate reads (mapped per gene) per sample from either HybPhaser or HybPiper folders
	if [[ $HYB == "phaser" ]]; then
		
		cat $SAMPLESDIR/$SAMPLE/reads/*_unpaired.fasta $SAMPLESDIR/$SAMPLE/reads/*_combined.fasta > $OUTDIR/$SAMPLE"_mapped_reads.fasta" 2>/dev/null
		
	else 
	
		cat $SAMPLESDIR/$SAMPLE/*/*_unpaired.fasta $SAMPLESDIR/$SAMPLE/*/*_interleaved.fasta > $OUTDIR/$SAMPLE"_mapped_reads.fasta" 2>/dev/null
	
	fi	
	
# The following pipline is based on Pierre Lindenbaum's script: https://www.biostars.org/p/3003/#3008 and composes of the following actions:
# in every row that starts with '>' put '@' at the end of the row
# replace '>' with '#'
# remove the break lines, replace '#' with a breakline, replace '@' with a tab
# sort the file according to the second column (sequences). the -u option keeps only one copy of each unique string.
# add '>' at the begining of each line
# sustitute tab (\t)  with a breakline (\n)
# remove the first line (it's a blank line) and save the result into $output

	if [[ $DEDUPSEQ == "TRUE" ]]; then 
		echo Deduplicating reads files for $SAMPLE

		sed -e '/^>/s/$/@/' -e 's/^>/#/' $OUTDIR/$SAMPLE"_mapped_reads.fasta" |\
		tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
		sort -u -t $'\t' -f -k 2,2  |\
		sed -e 's/^/>/' -e 's/\t/\n/' |\
		tail -n +2 > $OUTDIR/$SAMPLE"_mapped_reads_dedup-seq.fasta"

		rm $OUTDIR/$SAMPLE"_mapped_reads.fasta" 2>/dev/null
		
	else
		# NOTE: if the following commented code is run on paired-end interleaved data,
		#		which is the normal output from HybPhaser, one of the paired reads
		#		will be dropped (they have identical names), which seems like non-ideal
		#		behaviour (and different from what happens in removing duplicates above)

		#sed -e '/^>/s/$/@/' -e 's/^>/#/' $OUTDIR/$SAMPLE"_mapped_reads.fasta" |\
		#tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
		#sort -u -t $'\t' -f -k 1,1  |\
		#sed -e 's/^/>/' -e 's/\t/\n/' |\
		#tail -n +2 > $OUTDIR/$SAMPLE"_mapped_reads_nodupnames.fasta"
		
		#mv $OUTDIR/$SAMPLE"_mapped_reads_nodupnames.fasta" $OUTDIR/$SAMPLE"_mapped_reads.fasta"
	fi

done
