########################################################################
#   USAGE:
#       python3 StringSpelledByGenomePath.py
#   DESCRIPTION:
#       Reconstruct a string from consecutive k-mers of a genome path
#   BIOINFORMATICS II GENOME SEQUENCING:
#       1.3 String Reconstruction as a Walk in the Overlap Graph – Step 3
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20160525
########################################################################
import sys
# ---reconstruct_string-------------------------------------------------
# Reconstructs a string from consecutive k-mers of a genome path
# @param list of kmers
# @return reconstructed string
# ----------------------------------------------------------------------
def reconstruct_string(kmers):
# Initializes string to be returned with a reconstructed dna string
    reconstructedString = ''
    reconstructedString += kmers[0]
# Compares ith kmer's last k-1 symbols to ith+1 kmer's first k-1 symbols
    for i in range(len(kmers) - 1):
        if kmers[i][-(len(kmers[0]) - 1):] == kmers[i+1][:(len(kmers[0]) - 1)]:
            reconstructedString += kmers[i+1][-1:]
# Returns output
    return reconstructedString
# ----------------------------------------------------------------------
# ---MAINCODE-----------------------------------------------------------
kmers = sys.stdin.read().splitlines()
print (reconstruct_string(kmers))
# ----------------------------------------------------------------------
