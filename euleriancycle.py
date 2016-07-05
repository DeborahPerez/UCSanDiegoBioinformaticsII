########################################################################
#   USAGE:
#       python3 Euleriancycle.py
#   DESCRIPTION:
#       Construct an Eulerian cycle from an adjacency list of an
#       Eulerian directed graph
#   BIOINFORMATICS II GENOME SEQUENCING:
#       1.2 From Euler's Theorem to an Algorithm for Finding Eulerian
#       Cycles - Step 2
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20160619
########################################################################
import sys
# ---construct_eulerian_cycle-------------------------------------------
# Constructs an Eulerian cycle from an adjacency list of an Eulerian
# directed graph
# @param Eulerian directed graph
# @return Eulerian cycle
# ----------------------------------------------------------------------
def construct_eulerian_cycle(gMap):
    node = [next(iter(gMap))]
    cycle = []
    for i in node:
        cycle.append(i)
    print (node, cycle, gMap)
    while len(gMap) > 0:
        edges = gMap[node[0][0]]
        nextNode = edges[0]
        cycle.append(nextNode)
        if len(edges) == 1:
            del gMap[node[0]]
        else:
            del edges[0]

    return cycle
# ----------------------------------------------------------------------
# ---MAINCODE-----------------------------------------------------------
rawgMap = sys.stdin.read().splitlines()
gMap = {}
for i in range(len(rawgMap)):
    delimiter = ' -> '
    nodes = rawgMap[i].split(delimiter)
    for j in range(1, len(nodes)):
        delimiter = ','
        moreNodes = nodes[j].split(delimiter)
        gMap[nodes[0]] = moreNodes
print (construct_eulerian_cycle(gMap))
# ----------------------------------------------------------------------
# sample output
# 6->8->7->9->6->5->4->2->1->0->3->2->6
