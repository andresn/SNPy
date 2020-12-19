#!/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 4gb

###
# Find sequences aligned to a reference sequence with only SNPs (no indels)
# Note: run using snpy.sh
#
# E.g.,               to run: sbatch snpy.sh -r 'REF' -f test.fasta
#       to start at a column: sbatch snpy.sh -r 'REF' -f test.fasta -s 2
#                   for help: python3 snpy.py -h
#               to run tests: sbatch snpy.sh -t -r 'REF' -f test.fasta
###

module unload miniconda2
module load miniconda3
module load muscle

while getopts t:r:f: flag
do
    case "${flag}" in
        t) run_tests="-t";;
        r) ref_string="${OPTARG}";;
        f) fasta_file="${OPTARG}";;
    esac
done
python3 snpy.py $run_tests -r "$ref_string" -f "$fasta_file"