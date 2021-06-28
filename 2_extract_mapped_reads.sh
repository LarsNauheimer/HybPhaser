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
	
	SAMPLES='ls -d -1 $SAMPLESDIR/*'
	SAMPLE=${SAMPLES##*/}

fi


for i in $SAMPLESDIR/*
do
	SAMPLE=${i##*/}
	
	echo Generating deduplicated on-target reads files for $SAMPLE
	
	# concatenate reads (mapped per gene) per sample from either HybPhaser or HybPiper folders
	if [[ $HYB == "phaser" ]]; then
		
		cat $SAMPLESDIR/$SAMPLE/reads/*_unpaired.fasta $SAMPLESDIR/$SAMPLE/reads/*_combined.fasta > $OUTDIR/$SAMPLE"_mapped_reads_with_dups.fasta" 2>/dev/null
		
	else 
	
		cat $SAMPLESDIR/$SAMPLE/*/*_unpaired.fasta $SAMPLESDIR/$SAMPLE/*/*_interleaved.fasta > $OUTDIR/$SAMPLE"_mapped_reads_with_dups.fasta" 2>/dev/null
	
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
	
		sed '/^>/s/$/@/' < $OUTDIR/$SAMPLE"_mapped_reads_with_dups.fasta" | 
		sed 's/^>/#/' |\
		tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
		sort -u -f -k 2,2  |\
		sed -e 's/^/>/' |\
		tr '\t' '\n/' |\
		tail -n +2 > $OUTDIR/$SAMPLE"_mapped_reads_dedup-seq.fasta"	
		
	else
	
		sed '/^>/s/$/@/' < $OUTDIR/$SAMPLE"_mapped_reads_with_dups.fasta" | 
		sed 's/^>/#/' |\
		tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
		sort -u -f -k 1,1  |\
		sed -e 's/^/>/' |\
		tr '\t' '\n/' |\
		tail -n +2 > $OUTDIR/$SAMPLE"_mapped_reads.fasta"
		
	fi

	rm $OUTDIR/$SAMPLE"_mapped_reads_with_dups.fasta" 2>/dev/null
	
done
