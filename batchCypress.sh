#!/bin/bash
# batchCypress.sh
#
# usage:
# sh batchCypress.sh <input.list>
#
# Runs Cypress.sh
# in batch on genesis
# 

# Control Panel -----------------------

	# Input List
	INLIST=$1

	QSUB='qsub -q centos5.q -S /bin/bash -V -l  mem_token=8G'

	# Exon GTF File
	EXONS="$PWD/IRF5.gtf"

	# Genome Input Fasta
	GENOME="/home/ababaian/resources/genome/hg19r.fa"

	# This Directory
	HOMEDIR=$(pwd)
	
# Loop ================================

# Number of Iterations
	NUM=$(wc -l $INLIST | cut -f1 -d' ' - )

for X in $(seq $NUM)
do
	echo " ===== BATCH LOOP ===== "
	echo "       Iteration: $X"

	# Extract Line form input.list
	LINE=$(sed -n "$X"p $INLIST)

	# Output name
	NAME=$(echo $LINE | cut -f1 -d' ')

	# File Path
	FPATH=$(echo $LINE | cut -f2 -d' ')

	# Read Length
	READLEN=$(echo $LINE | cut -f3 -d' ')
	
	echo "        Sample: $NAME"
	echo "        File: $FPATH"
	echo "        Read Length: $READLEN"
	
	# Run Cypress Branch on node
	echo " Cufflinks Analysis"
	echo "    Cyp_branch.sh $FPATH $NAME $EXONS $GENOME $READLEN $HOMEDIR"
	
	# Genesis/Apollo
	$QSUB Cyp_branch.sh $FPATH $NAME $EXONS $GENOME $READLEN $HOMEDIR
	
	# Wait 2 minutes to offset genome reading
	sleep 120s

	# Local Submission
	#sh runCuff.sh $FPATH $NAME

	echo ""
	echo " ==================== "
	echo ""
	
done

#End of script
