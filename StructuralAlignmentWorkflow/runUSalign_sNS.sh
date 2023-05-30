#!/bin/bash
# To Run one to one pdb, outputting the translation and rotation matrix to a 
# output directory and file name
# usage: sh runUSalign_sNS.sh <query_struct> <target_struct> <out>
# 		query_struct - structure from library to be aligned to the 
#				target structure
#		target_struct- structure to be aligned to; of interest
#		out - path string used for naming output files
# USalign parameters:
#	-ter 2	(only consider the first chain of both structures to align to)
#	-mm 6 (use the semi-non-sequential alignment algorithm)
#	-outfmt 0 (provide output in full, human-readable form)
# 	-m $out.dat (write the translation and rotation matrix needed to 
#			recreate the alignment of query_struct to target_struct
#			to $out.dat)

# NOTE: avoid this line by adding the USalign executable within the python path/executable path when creating the environment
USALIGN_HOME=~/Apps/USalign

query_struct=$(readlink -f $1) 
target_struct=$(readlink -f $2)
out=$3

$USALIGN_HOME/USalign $query_struct $target_struct -ter 2 -mm 6 -outfmt 0 -m $out.dat
cat $out.dat
rm $out.dat

