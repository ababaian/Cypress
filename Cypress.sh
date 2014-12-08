#!/bin/bash
#
# Usage: Cypress.sh <Library.bam> <exons.gtf> <genome.fa> <read_length_opt.>
#
#

# Initialization ~~~~~~~~~~~~~~~~~~~~~~

LIBRARY=$1 # Library

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

READLEN=$4

if [ -z "$READLEN" ]
then
	echo " No read length set, using default 75"
	READLEN=75
fi


