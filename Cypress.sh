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
echo "         Library: $LIBRARY"
echo "         Exon File: $EXONS"
echo "         Genome: $GENOME"
echo "         Readlength: $READLEN"
echo "         Current Directory: $PWD"
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

# Samtools localization
# use local
	samtools='/home/ababaian/bin/samtools'

# Index Jucntions
# Splice Alignment -------------------------------------------------------------
# Align reads in Library to the splice junctions generated from Splice O Matic
#

# Generate the fastq file if it doesn't exist
if [ ! -f $NAME.1.fq ]
then
	echo ' Fastq file not found for the library so its goign to be generated'

	# Sort the library (Require Samtools 1.1)
	#$samtools sort -n -O 'bam' -T tempSort $LIBRARY > tempSort.bam
	
	# Samtools pre 1.1
	samtools sort -n $LIBRARY tempSort

	echo ' Post sort ls'
	ls -alh
	
	if [ ! -f tempSort.bam ]
	then
		echo ' tempSort.bam not generated! Error 3!'
		exit 3
	fi

	# Convert BAM file to FASTQ
	bam2fastx -q -Q -A -P -N -o $NAME.fq tempSort.bam
fi

# Align reads to the the junctions
	bowtie2-build junction.fa junction

	# High fidelity setting
	bowtie2 -x junction -p 2 -1 $NAME.1.fq -2 $NAME.2.fq --end-to-end -a -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --no-unal -S $NAME.jnc.sam

	# Normal Fidelity
	#bowtie2 -x junction -p 2 -1 $NAME.1.fq -2 $NAME.2.fq --end-to-end -a -D 15 -R 2 -L 22 -i S,1,1.15 --no-unal -S $NAME.jnc.sam

	# Normal Fidelity only primary alignments
	#bowtie2 -x junction -p 2 -1 $NAME.1.fq -2 $NAME.2.fq --end-to-end -a -D 15 -R 2 -L 22 -i S,1,1.15 --no-unal -S $NAME.jnc.sam


	echo ' Post alignment ls'
	ls -alh

	# Convert sam file into a bam file
	$samtools sort -O 'bam' -T tempspace.tmp -o junction.bam $NAME.jnc.sam

	# Retain only primary reads
	$samtools view -F 256 junction.bam -b -o junction.primary.bam
	$samtools index junction.primary.bam

	echo ' Post processing/sorting ls'
	ls -alh
# Analysis of Alignments -------------------------------------------------------

cut -f4 junction.bed | awk -v RL="$READLEN" '{ \
	print $1 "\t"\
	RL-5 "\t"\
	RL+5 "\t"\
	$1 "\t"\
	"0\t"\
	"+\t"\
 	"0\t"\
	"0\t"\
	"255,0,0\t"\
	"1\t"\
	"10\t"\
	"0\t"}' - > jncCenters.bed.tmp

bedtools coverage -abam junction.primary.bam -b jncCenters.bed.tmp | cut -f1,13 - | sed 's/_/\t/g' - > junctionCoverage.tmp

# awk -v readlen="$READLEN" '{print $1 "\t" readlen - 5 "\t" readlen + 4  "\t" $1}' - 

# Splice Donors
cut -f1 junctionCoverage.tmp | sort -u - > doner.tmp
dLen=$(wc -l doner.tmp | cut -f1 -d' ' - ) 

# Splice Acceptors
cut -f2 junctionCoverage.tmp | sort -u - > accept.tmp
aLen=$(wc -l accept.tmp | cut -f1 -d' ' - )

# Initilize output matrix
	rm -f jncMatrix.txt
	rm -f jncMatrix.tmp
	touch jncMatrix.tmp

# Iterate through each splice doner as a new line
for Nd in $(seq 1 $dLen)
do

	# Initilize output Row for this doner
	rm -f jncRow.tmp
	touch jncRow.tmp

	# Extract Splice Doner Name
	DonerName=$(cut -f1 doner.tmp | sed -n "$Nd"p -)

	# Iterate through each splice accepter as new column
	for Na in $(seq 1 $aLen)
	do
		# Extract Splice Acceptor Name
		AcceptName=$(cut -f1 accept.tmp | sed -n "$Na"p -)

		# Extract the number of reads over the splice junction
		LookUp=$(echo "$DonerName\t$AcceptName")
		SpliceScore=$(grep -P "$LookUp" junctionCoverage.tmp | cut -f3 - )

		if [ "$SpliceScore" = '' ]
		then
			SpliceScore='0'
		fi

		# Append to row
		cp jncRow.tmp holder.tmp
		echo $SpliceScore | paste holder.tmp - > jncRow.tmp

		#echo $DonerName $AcceptName $SpliceScore
		#cat jncRow.tmp
		
	done
	cat jncRow.tmp

	# Append completed row to matrix output
	cat jncRow.tmp >> jncMatrix.txt
done

mv doner.tmp  doner.txt
mv accept.tmp accept.txt 
rm *.tmp

# End of script ; )
