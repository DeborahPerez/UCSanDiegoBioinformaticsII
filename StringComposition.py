########################################################################
#   USAGE:
#       python3 StringComposition.py
#   DESCRIPTION:
#       Find all possible k-mers including repeats in a string.
#   BIOINFORMATICS II GENOME SEQUENCING:
#       1.2 The String Reconstruction Problem – Step 3
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20160525
########################################################################
import sys
# ---string_composition-------------------------------------------------
# Finds all possible k-mers including repeats in a string.
# @param dna and kmer length 'k'
# @return list of kmers on separate lines
# ----------------------------------------------------------------------
def string_composition(kValue, Dna):
# Initializes list to be returned with all possible kmers and repeats
    allKmers = []
# Slides window across the dna string while appending kmers to a list
    for i in range(len(Dna) - kValue + 1):
        allKmers.append(Dna[i:i+kValue])
# Return output
    return allKmers
# ----------------------------------------------------------------------
# ---MAINCODE-----------------------------------------------------------
kValueAndDna = sys.stdin.read().splitlines()
kValue = int(kValueAndDna[0])
Dna = kValueAndDna[1]
print ("\n".join(string_composition(kValue, Dna)))
# ----------------------------------------------------------------------
