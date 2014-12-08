# splice-O-matic 5000
# 
# Usage: spliceOmatic.sh <exons.gtf> <genome.fa> <read_length (opt.)>
#
# Output: <junctions.bed> <junctions.fa>
#
# The splice-O-matic offers the cutting edge in generating splice junctions
# between exons. What it doesn't do elegantly it does through pure brute
# mechanical force.
#
# For every reasonable exon pair in exon.bed the SOM will generate splice
# junctions and append them to a splice junction fasta file 
# 
# This can then be used to do a detailed search of an RNAseq file
# to quantify the relative abundance of each splice junction in a library
#

# File Format for Exons.gtf 
#chr	annotation	exon_position	start	end	exon_num	strand	NA	name
#		   (5end / int_exon / 3end)

# Splice Rules =================================================================
#===============================================================================
# 1) Exon 1 do not have an upstream splice pair, only a downstream one
# 2) Last Exons (3end) do not have a downstream splice pair, only an upstream one
# 3) All other splice combinations will be considered.
#===============================================================================

# Initialization ~~~~~~~~~~~~~~~~~~~~~~

EXONS=$1 # Exon File

if [ ! -f $EXONS ]
then
	echo " Warning, the <exons.gtf> file $EXONS, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

GENOME=$2 # Genome File

if [ ! -f $GENOME ]
then 
	echo " Warning, the <genome.fa> file $GENOME, was not found"
	echo " Malfunction imminant"
	echo " shutting down!"
	exit 1
fi

READLEN=$3

if [ -z "$READLEN" ]
then
	echo " No read length set, using default 75"
	READLEN=75
fi

echo " Splice-O-Matic 5000 is awakening! Such glory!"
echo "     Run Parameters:"
echo "         Exon File: $EXONS"
echo "         Genome: $GENOME"
echo "         Readlength: $READLEN"
echo ''

# Script Core ==================================================================

# Partiion exons.gtf
#	Remove commented lines in exons GTF file if they exist
#	change the GTF score column to the exon-length
#	Flag if exon is larger then readlength (1) or not (0) 
# 	split the file into 5ends, internal exons, 3ends

	# Parse + calculate exon size
	sed '/^#/d' $EXONS | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5 - $4"\t"$7"\t"$8"\t"$9}' - > exon.tmp

	# Flag small exons
	awk -v readlen="$READLEN" '{if ($6 >= readlen) print 1;  else  print 0 }' exon.tmp > flag.tmp
	cut -f1-7 exon.tmp > start.tmp
	cut -f9 exon.tmp | paste start.tmp flag.tmp - > exon.final.tmp

	rm start.tmp flag.tmp exon.tmp 

	# Extract 3 ends
	grep '5end' exon.final.tmp > 5end.tmp
	# Extract internal exons
	grep 'int' exon.final.tmp > int.tmp
	# Extract 3 ends
	grep '3end' exon.final.tmp > 3end.tmp

	# Which chromosome is this gene on?
	CHR=$(sed '/^#/d' $EXONS | head -n1 - | cut -f1 - )
	# Which strand is the gene on?
	STRAND=$(sed '/^#/d' $EXONS | head -n1 - | cut -f7 - )

# Build Splice Junction Map
#
#	0) Check all exons are greater then read-length
#	1) For every 5end exon, connect to each downstream exon
#	2) For every intexon, connect to each downstream exon
#	3) Format a master splice file

# 0) Check exon lengths
# for early alpha, treat small exons the same as larger exons but output an error

# 1) 5end Exon Builds
#

# Loop through exon files
	ITER_5=$(wc -l 5end.tmp | cut -f1 -d' ' -)
	ITER_i=$(wc -l int.tmp | cut -f1 -d' ' -)
	ITER_3=$(wc -l 3end.tmp | cut -f1 -d' ' -)

# Splice junctions from the 5' end
for N in $(seq 1 $ITER_5)
do
	# Extract the Nth entry in 5end.tmp for splice doner
	LINE=$(sed -n $N'p' 5end.tmp)
	
	#Splice Doner
	DONER=$(echo $LINE | cut -f5 -d' ' -)
	DNAME=$(echo $LINE | cut -f9 -d' ' -)

	# Loop through int.tmp and 3end.tmp for splice acceptors
	#Int loop
	for M in $(seq 1 $ITER_i)
	do
		#Extract Nth Acceptor site in each entry in int.tmp
		ACCEPT=$(sed -n $M'p' int.tmp | cut -f4 - )
		ANAME=$(sed -n $M'p' int.tmp | cut -f9 - )

		if [ $DONER -lt $ACCEPT ]
		then
		# Junction possible
			# Append Doner~Acceptor to output splice file
			echo $DONER\~$ACCEPT\;$DNAME\_$ANAME >> splice_junction.tmp
			
		#else
		# Junction not possible skip
		fi
	done

	#3end loop
	for M in $(seq 1 $ITER_3)
	do
		#Extract Nth Acceptor site in each entry in int.tmp
		ACCEPT=$(sed -n $M'p' 3end.tmp | cut -f4 - )
		ANAME=$(sed -n $M'p' 3end.tmp | cut -f9 - )

		if [ $DONER -lt $ACCEPT ]
		then
		# Junction possible
			# Append Doner~Acceptor to output splice file
			echo $DONER\~$ACCEPT\;$DNAME\_$ANAME >> splice_junction.tmp
		#else
		# Junction not possible skip
		fi
	done
done

# Splice junctions from the internal exons
for N in $(seq 1 $ITER_i)
do
	# Extract the Nth entry in int.tmp for splice doner
	LINE=$(sed -n $N'p' int.tmp)
	
	#Splice Doner
	DONER=$(echo $LINE | cut -f5 -d' ' -)
	DNAME=$(echo $LINE | cut -f9 -d' ' -)

	# Loop through int.tmp and 3end.tmp for splice acceptors
	#Int loop
	for M in $(seq 1 $ITER_i)
	do
		#Extract Nth Acceptor site in each entry in int.tmp
		ACCEPT=$(sed -n $M'p' int.tmp | cut -f4 )
		ANAME=$(sed -n $M'p' int.tmp | cut -f9 - )

		if [ $DONER -lt $ACCEPT ]
		then
		# Junction possible
			# Append Doner~Acceptor to output splice file
			echo $DONER\~$ACCEPT\;$DNAME\_$ANAME >> splice_junction.tmp

		#else
		# Junction not possible skip
		fi
	done

	#3end loop
	for M in $(seq 1 $ITER_3)
	do
		#Extract Nth Acceptor site in each entry in int.tmp
		ACCEPT=$(sed -n $M'p' 3end.tmp | cut -f4 - )
		ANAME=$(sed -n $M'p' 3end.tmp | cut -f9 - )

		if [ $DONER -lt $ACCEPT ]
		then
		# Junction possible
			# Append Doner~Acceptor to output splice file
			echo $DONER\~$ACCEPT\;$DNAME\_$ANAME >> splice_junction.tmp
		#else
		# Junction not possible skip
		fi
	done
done

# Remove Duplicate Entries
# Format for awk
	sort -u -t\; -k1,1 splice_junction.tmp > jnc.u.tmp
	sed 's/~/\t/g' jnc.u.tmp | sed 's/\;/\t/g' - > jnc.tmp

# Convert Junction coordinates (1 base) to BED12 foramt (0 base)

awk -v RL="$READLEN" -v CHR="$CHR" -v STRAND="$STRAND" '{ \
	print CHR "\t"\
	$1-RL "\t"\
	$2+RL "\t"\
	$3 "\t"\
	"0\t"\
	STRAND"\t"\
	"0\t"\
	"0\t"\
	"255,0,0\t"\
	"2\t"\
	RL","RL"\t"\
	"0," ($2-$1+RL) "\t"}' jnc.tmp > junction.bed

# Extract sequences from genome (Master splice file)
	fastaFromBed -split -fi $GENOME -bed junction.bed -fo junction.fa.tmp

# Interlace the exon junction names into the fasta output file

	# Fasta Sequences (even numbers)
	sed '1d; n; d' junction.fa.tmp > seq.fa.tmp

	# Junction Names
	cut -f3 jnc.tmp | sed 's/^/>/g' - | paste -d'\n' - seq.fa.tmp > junction.fa

#Cleanup
	rm *.tmp




