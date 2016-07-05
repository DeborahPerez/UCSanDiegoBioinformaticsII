########################################################################
#   USAGE:
#       python3 DeBruijnGraphFromkmers.py
#   DESCRIPTION:
#       Construct the de Bruijn graph from a set of k-mers
#   BIOINFORMATICS II GENOME SEQUENCING:
#       1.5 Walking in the de Bruijn Graph – Step 7
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20160605
########################################################################
# ---kmer_debruijn_graph------------------------------------------------
# Constructs a de Bruijn graph from a set of k-mers
# @param list of kmers
# @return deBruijn graph
# ----------------------------------------------------------------------
import sys
def kmer_debruijn_graph(kmers):
    graphMap = {}
    for kmer in kmers:
        if _prefix(kmer) not in graphMap:
            graphMap[_prefix(kmer)] = [_suffix(kmer)]
        else:
            graphMap[_prefix(kmer)].append(_suffix(kmer))
    return graphMap
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
kmers = sys.stdin.read().splitlines()
graphMap = kmer_debruijn_graph(kmers)
for key in sorted(graphMap):
    print (key + ' -> ' + ",".join(graphMap[key]))
# ----------------------------------------------------------------------
