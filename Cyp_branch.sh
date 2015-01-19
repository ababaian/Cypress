#!/bin/bash
# Cypress call-script for batch runs
# Usage: Cyp_branch.sh <Library.bam> <exons.gtf> <genome.fa> <name> <read_length_opt.> <parent_dir>
#
#

# Initialization ~~~~~~~~~~~~~~~~~~~~~~

LIBRARY=$1 # Library

NAME=$2

if [ ! -f $LIBRARY ]
then
	echo " Warning, the <Library.bam> file $LIBRARY, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

EXONS=$3 # Exon File

if [ ! -f $EXONS ]
then
	echo " Warning, the <exons.gtf> file $EXONS, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

GENOME=$4 # Genome File

if [ ! -f $GENOME ]
then 
	echo " Warning, the <genome.fa> file $GENOME, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

READLEN=$5 # Read Length

if [ -z "$READLEN" ]
then
	echo " No read length set, using default 75"
	READLEN=75
fi

# Parent Directory
parDIR=$6

echo " Cypress initializing with parameters"
echo "         Name: $NAME"
echo "         Exon File: $EXONS"
echo "         Genome: $GENOME"
echo "         Readlength: $READLEN"
echo "         Parent Dir: $parDIR"
echo ''

# Script Core ==================================================================
# This script is executed on an execution node of a computer
# Initilize it to run Cypress.sh with what it needs!

# Make a temporary directory
	mkdir /tmp/Cypress/

# Go to parent directory
	cd $parDIR

# Copy core files for running Cypress to the node
	cp Cypress.sh /tmp/Cypress/Cypress.sh
	cp spliceOmatic.sh /tmp/Cypress/spliceOmatic.sh
	cp -L $EXONS /tmp/Cypress/input.gtf
	cp -L $LIBRARY /tmp/Cypress/$NAME.bam

# Move to temporary directory
	cd /tmp/Cypress
	samtools index $NAME.bam

	echo " Cypress Branch has been initialized on the node"
	echo " Folder contents are:"
	ls -alh
	echo ''
	echo ''


# Run Cypress locally (finaly)
echo "Running Cypress:"
echo "  sh Cypress.sh $NAME.bam input.gtf $GENOME $READLEN"
sh Cypress.sh $NAME.bam input.gtf $GENOME $READLEN 

echo " Post Cypress LS"

ls -alh

# Initilize output directory
	mkdir $parDIR/$NAME

# Move ouput 
	cp junction.fa $parDIR/$NAME/
	cp junction.bam $parDIR/$NAME/
	cp junction.bam.bai $parDIR/$NAME/
	cp junction.bed $parDIR/$NAME/
	cp junction.primary.bam $parDIR/$NAME/
	cp junction.primary.bam.bai $parDIR/$NAME/
	cp jncMatrix.txt $parDIR/$NAME/
	cp doner.txt $parDIR/$NAME/
	cp accept.txt $parDIR/$NAME/

echo " batch job of Cypress completed"
