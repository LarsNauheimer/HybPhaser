#!/bin/bash
############
set -e
set -u
set -o pipefail

NAMESLIST="nameslist.txt"
BASE_FOLDER="./"
OUTPUT_FOLDER=$BASE_FOLDER/'_mapped_reads'

while getopts 'n:b:o:' OPTION; do
  case "$OPTION" in

    n)
       NAMESLIST=$OPTARG
     ;;
    
    b)
      BASE_FOLDER=$OPTARG
     ;;         
    
    o)
      OUTPUT_FOLDER=$OPTARG
     ;;         
    
    ?)
      echo "script usage: $(basename $0) [-n file with samples names (nameslist.txt)] [-b HybPiper base folder] [-o output folder]"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

mkdir -p $OUTPUT_FOLDER/

while read NAME;
do
if [ $NAME != "" ]; then 
    echo Extracting mapped reads from $NAME
    #samtools fastq -F 4 $BASE_FOLDER/$NAME/$NAME".bam" > $OUTPUT_FOLDER/$NAME".fastq"
    samtools view -F 4 -b $BASE_FOLDER/$NAME/$NAME".bam"| samtools sort -n -o $OUTPUT_FOLDER/$NAME"_sort.bam"
    samtools fixmate -m $OUTPUT_FOLDER/$NAME"_sort.bam" $OUTPUT_FOLDER/$NAME"_fixm.bam"
    samtools sort -o $OUTPUT_FOLDER/$NAME"_fixmsort.bam" $OUTPUT_FOLDER/$NAME"_fixm.bam"
    samtools markdup -r -s $OUTPUT_FOLDER/$NAME"_fixmsort.bam" $OUTPUT_FOLDER/$NAME"_dupm.bam"
    samtools fastq $OUTPUT_FOLDER/$NAME"_dupm.bam" > $OUTPUT_FOLDER/$NAME".fastq"
    rm $OUTPUT_FOLDER/$NAME"_sort.bam" $OUTPUT_FOLDER/$NAME"_fixm.bam" $OUTPUT_FOLDER/$NAME"_fixmsort.bam" $OUTPUT_FOLDER/$NAME"_dupm.bam"
fi
done < $NAMESLIST




