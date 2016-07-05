########################################################################
#   USAGE:
#       python3 DeBruijnGraph.py
#   DESCRIPTION:
#        Solve the De Bruijn Graph from a String Problem
#   BIOINFORMATICS II GENOME SEQUENCING:
#       1.4 Another Graph for String Reconstruction – Step 6
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20160604
########################################################################
import sys
# ---construct_de_bruijn_graph------------------------------------------
# Reconstructs a string from consecutive k-mers of a genome path
# @param an integer k and a string Text.
# @return  DeBruijn(Text), in the form of an adjacency list.
# ----------------------------------------------------------------------
def construct_de_bruijn_graph(kValue, dna):
    kMinus1mers = string_composition(kValue-1, dna)
    graphMap = {}
    for p in range(len(dna) - kValue + 1):
        if kMinus1mers[p] not in graphMap:
            graphMap[kMinus1mers[p]] = [kMinus1mers[p+1]]
        else:
            graphMap[kMinus1mers[p]].append(kMinus1mers[p+1])
    return graphMap
# ---string_composition-------------------------------------------------
# Finds all possible k-mers including repeats in a string.
# @param dna and kmer length 'k'
# @return list of kmers on separate lines
# ----------------------------------------------------------------------
def string_composition(kValue, dna):
# Initializes list to be returned with all possible kmers and repeats
    allKmers = []
# Slides window across the dna string while appending kmers to a list
    for i in range(len(dna) - kValue + 1):
        allKmers.append(dna[i:i+kValue])
# Return output
    return allKmers
# ----------------------------------------------------------------------
# ---_prefix------------------------------------------------------------
# @param text
# @return prefix
# ----------------------------------------------------------------------
def _prefix(text):
    return text[:len(text)-1]
# ----------------------------------------------------------------------
# ---_suffix------------------------------------------------------------
# @param text
# @return suffix
# ----------------------------------------------------------------------
def _suffix(text):
    return text[-(len(text) - 1):]
# ----------------------------------------------------------------------
# ---MAINCODE-----------------------------------------------------------
kValueAndDna = sys.stdin.read().splitlines()
kValue = int(kValueAndDna[0])
dna = kValueAndDna[1]
graphMap = construct_de_bruijn_graph(kValue, dna)
for key in sorted(graphMap):
    print (key + ' -> ' + ",".join(graphMap[key]))
# ----------------------------------------------------------------------
