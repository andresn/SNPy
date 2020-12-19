###
# Find sequences aligned to a reference sequence with only SNPs (no indels)
# Note: run using snpy.sh
###

import argparse
from time import time 
import json
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from io import StringIO

parser = argparse.ArgumentParser()
parser.add_argument('-t', action='store_true', help='Runs unit tests based off of ./test.aln.')
parser.add_argument('-r', action='store', dest='ref_string', type=str, help='The entire string after ">" in your FASTA file describing your reference sequence.')
parser.add_argument('-f', action='store', dest='fasta_file', type=str, help='The .fasta file containing the sequences you want aligned and analyzed.')
parser.add_argument('-s', action='store', dest='start_index', type=str, help='The index of the column where you want to start the analysis from.')

args = parser.parse_args()
run_tests = args.t
ref_string = args.ref_string.split(' ')[0]
fasta_file = args.fasta_file
start_index = int(args.start_index) if args.start_index != '' else 0

time_stamp = str(time()).split('.')[0][-7:]

if run_tests:
    print('Running tests...')
    alignment = AlignIO.read('test.aln', 'clustal')
else:
    out_file = 'snpy-' + time_stamp + '.aln'
    cline = MuscleCommandline('muscle', input=fasta_file, clw=True, out=out_file)
    cline()
    alignment = AlignIO.read(out_file, 'clustal')

indel_hash = {}
snp_hash = {}

total_columns = alignment.get_alignment_length()

column_index = 3
row_index = 1
total_rows = len(alignment)

print('ref_string: ', ref_string)

for i in range(0, total_rows):
    print(alignment[i].id)
    if alignment[i].id == ref_string:
        ref_row_index = i

def parse_alignment_column(alignment, indel_hash, snp_hash, total_columns, total_rows, start_index=0):
    for column_index in range(start_index, total_columns):
        for row_index in range(0, total_rows):
            if alignment[row_index, column_index] != alignment[0, column_index]:
                if alignment[row_index].id not in indel_hash:
                    if alignment[0, column_index] == "-" or alignment[row_index, column_index] == "-":
                        indel_hash[alignment[row_index].id] = True
                    else:
                        snp_hash[alignment[row_index].id] = True
    return indel_hash, snp_hash

column_index = 0
indel_hash, snp_hash = parse_alignment_column(alignment, indel_hash, snp_hash, total_columns, total_rows, start_index)

indel_array = []
for hash, value in indel_hash.items() :
    if value:
        indel_array.append(hash)
SNP_report = []
for hash, value in snp_hash.items() :
    if value and hash not in indel_array:
        SNP_report.append(hash)

print('snpy.py job complete: see snpy-' + time_stamp + '.out and .aln')
if run_tests:
    print('SNPs with no indels in sequences: ' + (', ').join(SNP_report))
    print('indel_hash:')
    print(indel_hash)
    print('snp_hash:')
    print(snp_hash)
    assert indel_hash == {'825225': True, '833242': True}
    assert snp_hash == {'825225': True, '736945': True}
else:
    with open('snpy-' + time_stamp + '.out', 'w') as report_file:
        report_file.write(
            '\n'.join([
                'SNPs with no indels in sequences: ' + (', ').join(SNP_report),
                'indel_hash:',
                json.dumps(indel_hash),
                'snp_hash:',
                json.dumps(snp_hash)
            ])
        )
