# SNPy
Find sequences aligned to a reference sequence with only SNPs (no indels).

Note: run using ./snpy.sh, or sbatch snpy.sh etc. if on HPPC cloud, and "REF" is everything after ">" in the FASTA file where your reference sequence is:
```
E.g.:               to run: sbatch snpy.sh -r 'REF' -f test.fasta
      to start at a column: sbatch snpy.sh -r 'REF' -f test.fasta -s 2
                  for help: python3 snpy.py -h
              to run tests: sbatch snpy.sh -t -r 'REF' -f test.fasta
```
