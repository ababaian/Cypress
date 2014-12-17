#!/bin/bash
#
# Usage: Cypress.sh <Library.bam> <exons.gtf> <genome.fa> <read_length_opt.>
#
#

# Initialization ~~~~~~~~~~~~~~~~~~~~~~

LIBRARY=$1 # Library

NAME=$(echo $LIBRARY | sed 's/\.bam//g' - )

if [ ! -f $LIBRARY ]
then
	echo " Warning, the <Library.bam> file $LIBRARY, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

EXONS=$2 # Exon File

if [ ! -f $EXONS ]
then
	echo " Warning, the <exons.gtf> file $EXONS, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

GENOME=$3 # Genome File

if [ ! -f $GENOME ]
then 
	echo " Warning, the <genome.fa> file $GENOME, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

READLEN=$4 # Read Length

if [ -z "$READLEN" ]
then
	echo " No read length set, using default 75"
	READLEN=75
fi

echo " Cypress Analysis parameters valid."
echo " Running with parameters..."
echo "         Name: $NAME"
echo "         Exon File: $EXONS"
echo "         Genome: $GENOME"
echo "         Readlength: $READLEN"
echo ''

# Script Core ==================================================================

# Splice Junctions -------------------------------------------------------------
# Run Splice O Matic to generate junction.bed // junction.fa
if [ ! -f junction.bed ] | [ ! -f junction.fa ]
then
	echo " Running Splice O Matic ..."
	echo "      sh spliceOmatic.sh $EXONS $GENOME $READLEN "
	echo ''

	sh spliceOmatic.sh $EXONS $GENOME $READLEN
fi

if [ ! -f junction.bed ] || [ ! -f junction.fa ]
then 
	echo " Warning, either junction.bed or junction.fa was not generated"
	echo " Check for errors in SpliceOMatic"
	echo " Exiting"
	exit 2
fi

# Index Jucntions
# Splice Alignment -------------------------------------------------------------
# Align reads in Library to the splice junctions generated from Splice O Matic
#

# Generate the fastq file if it doesn't exist
if [ ! -f tempSort.1.fq ]
then
	# Sort the library (Require Samtools 1.1)
	samtools sort -n -O 'bam' -T tempSort $LIBRARY > tempSort.bam

	# Convert BAM file to FASTQ
	bam2fastx -q -Q -A -P -N -o $NAME.fq tempSort.bam
fi

# Align reads to the the junctions
	bowtie2-build junction.fa junction

	bowtie2 -x junction -p 2 -1 $NAME.1.fq -2 $NAME.2.fq --end-to-end -a -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --no-unal -S $NAME.jnc.sam
	#tophat2 -p 1 -r 50 --report-secondary-alignments -g 1000 -o $PWD $GENOME temp.1.1fq temp.2.fq









