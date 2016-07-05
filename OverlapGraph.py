########################################################################
#   USAGE:
#       python3 OverlapGraph.py
#   DESCRIPTION:
#       Construct a graph of overlappting k-mers
#   BIOINFORMATICS II GENOME SEQUENCING:
#       1.3 String Reconstruction as a Walk in the Overlap Graph
#       – Step 9
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20160530
########################################################################
import sys
# ---construct_overlap_graph--------------------------------------------
# Constructs a graph of overlapping kmers
# @param a collection Patterns of k-mers
# @return the overlap graph Overlap(Patterns), in the form of an
# adjacency list.
# ----------------------------------------------------------------------
def construct_overlap_graph(kmers):
    graph = []
    for i in range(len(kmers)):
        nodes = []
        for j in range(len(kmers)):
            if _suffix(kmers[i]) == _prefix(kmers[j]) and i != j:
                nodes.append(kmers[j])
        if len(nodes) != 0:
            nodeDescriptor = kmers[i] + ' -> ' + ",".join(nodes)
            graph.append(nodeDescriptor)
    return graph
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
# ---MAINCODE-----------------------------------------------------------
kmers = sys.stdin.read().splitlines()
overlapGraph = (construct_overlap_graph(kmers))
print ("\n".join(overlapGraph))
# ----------------------------------------------------------------------
