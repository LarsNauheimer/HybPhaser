#!/bin/bash
############


set -e
set -u
set -o pipefail

CONTIG="normal"
THREADS=1
SAMPLE=""
ALLELE_FREQ=0.15
READ_DEPTH=10
ALLELE_COUNT=4
BASE_FOLDER="./"

while getopts 'st:n:a:f:d:b:' OPTION; do
  case "$OPTION" in

    s)
      echo "supercontig"
      CONTIG="SC"
      ;;

    t)
      THREADS=$OPTARG
     ;;
    
    n)
    echo "name of sample=$OPTARG"
      SAMPLE=$OPTARG
     ;;
	 
    f)
      ALLELE_FREQ=$OPTARG
     ;;
    
    a)
      ALLELE_COUNT=$OPTARG
     ;;
    
    d)
      READ_DEPTH=$OPTARG
     ;;
    
    b)
      BASE_FOLDER=$OPTARG
     ;;         
    
    ?)
      echo "script usage: $(basename $0) [-s] supercontig [-t no of threads '1'] [-n name of sample] [-f  minimum allele frequency to use ambiguity code '0.1'] [-a min count of alleles to be regarded] [-d min coverage at site]" >&2
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

cd $BASE_FOLDER$SAMPLE
  for DIR in */
    do
    if [ -d "$DIR$SAMPLE" ]; then
            
      GENE=${DIR/\//}
      OUTDIR=$DIR$SAMPLE/remapping/
      OUTDIRCONSENSUS=$DIR$SAMPLE/sequences/remapping/

      # checking for paired-end reads
      if [ -f $DIR$GENE"_interleaved.fasta" ];then
        READS=$DIR$GENE"_interleaved.fasta"
        READTYPE="PE"
      else 
         READS=$DIR$GENE"_unpaired.fasta"
         READTYPE="SE"
      fi 
      
      if [ $CONTIG = 'SC'  ];then
       SC="_supercontig"
       REFFILE=$DIR$SAMPLE/sequences/intron/$GENE"_supercontig.fasta"
       OUTFILECONSENSUS=$OUTDIRCONSENSUS$GENE$"_supercontig-consensus.fasta"
       REF=$OUTDIR$GENE$SC".fasta"
      else 
       SC=""
       REFFILE=$DIR$SAMPLE/sequences/FNA/$GENE".FNA"
       OUTFILECONSENSUS=$OUTDIRCONSENSUS$GENE$"_consensus.fasta"
       REF=$OUTDIR$GENE".FNA"
      fi 


      BAM=$OUTDIR$GENE$SC"_remapped.bam"
      VCF=$OUTDIR$GENE$SC"_remapped.vcf"
      VCFZ=$OUTDIR$GENE$SC"_remapped.vcf.gz"

      mkdir -p $OUTDIR
      mkdir -p $OUTDIRCONSENSUS
      if [ -f "$REFFILE" ]; then
      
        cp "$REFFILE" "$OUTDIR"
      
        bwa index $REF
      
        if [ "$READTYPE" = "SE" ];then 
          bwa mem $REF $READS -t $THREADS | samtools sort > $BAM
        else 
          bwa mem -p $REF $READS -t $THREADS | samtools sort > $BAM
        fi
      
        bcftools mpileup -I -Ov -f $REF $BAM | bcftools call -mv -A -Oz -o $VCFZ
		bcftools index -f $VCFZ
        bcftools consensus -I -i "(DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= $ALLELE_FREQ && (DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= $READ_DEPTH && (DP4[2]+DP4[3]) >= $ALLELE_COUNT " -f $REF $VCFZ | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' > $OUTFILECONSENSUS
        echo "" >> $OUTFILECONSENSUS
      
        if [ "$CONTIG" != "s" ];then 
          sed -i "s/>.*/>$SAMPLE-$GENE/" $OUTFILECONSENSUS
        fi
      else 
        echo $REFFILE does not exist!
      fi
      
    else 
    echo $DIR$SAMPLE does not exist
    fi
    done
 cd ..

