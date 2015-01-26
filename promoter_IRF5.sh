#!/bin/bash
#
# Promoter_IRF5.sh <input jncMatrix.txt> <Library_name>
# Cypress Analysis II
# Promoter Contribution for IRF5

# INPUT
#INPUT=$1
INPUT='jncMatrix.txt'

# Library
#NAME=$2
NAME='A0123'

# IRF5 Native Promoter
# Rows: 1,3,4,6-9,23,24
# Cols: 2-6

# IRF5 Chimeric Promoter
# Rows: 5,21,22
# Cols: 2-6

# Generate Native Promoter Matrix
cut -f 2-5,14 $INPUT | sed -n 1p - > IRF5.native.tmp
cut -f 2-5,14 $INPUT | sed -n 3,4p - >> IRF5.native.tmp
cut -f 2-5,14 $INPUT | sed -n 6,9p - >> IRF5.native.tmp
cut -f 2-5,14 $INPUT | sed -n 23,24p - >> IRF5.native.tmp

# Calculate Native Matrix Sum
# echo $(RowSum) | ColSum
NatSum=$(echo $(sed 's/\t/+/g' IRF5.native.tmp | bc)  | sed 's/ /+/g' - | bc)

# Generate Chimeric Promoter Matrix
cut -f 2-5,14 $INPUT | sed -n 5p - > IRF5.chim.tmp
cut -f 2-5,14 $INPUT | sed -n 21,22p - >> IRF5.chim.tmp

ChiSum=$(echo $(sed 's/\t/+/g' IRF5.chim.tmp | bc)  | sed 's/ /+/g' - | bc)

echo $NatSum $ChiSum > output.txt
