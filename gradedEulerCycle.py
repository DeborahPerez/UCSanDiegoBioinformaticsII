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
#   VERSION:    20160730
########################################################################
import sys
# ---construct_eulerian_cycle-------------------------------------------
# Constructs an Eulerian cycle from an adjacency list of an Eulerian
# directed graph
# @param Eulerian directed graph
# @return Eulerian cycle
# ----------------------------------------------------------------------
def construct_eulerian_cycle(gMap):
    node = next(iter(gMap))
    cycle = []
    cycle.append(node)

# while loop for gMap with existing contents
    while len(gMap) > 0:
        if node in gMap:
            steps = gMap[node]
            step = steps[0]
            cycle.append(step)
#dictionary deletion mechanism
            if len(steps) > 1:
                del gMap[node][0]
            else:
                del gMap[node]
            node = step
#In the case that gMap is not empty and node is not in gMap, we will traverse
#the cycle and choose an index existing in gMap and continue from there
        else:
            for i in range(len(cycle) - 1):
                startPoint = cycle[i]
                if startPoint in gMap:
                    cutOffPoint = i
                    node = startPoint
#record last edge to compare to current cycle
            lastEdge = cycle[len(cycle) - 1]

#Modification and rearragement of cycle
            newCycle = []
            for edge in cycle[cutOffPoint:]:
                newCycle.append(edge)
            for moreEdge in cycle[1:cutOffPoint + 1]:
                newCycle.append(moreEdge)
            cycle = newCycle

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
cycle = construct_eulerian_cycle(gMap)
print ("->".join(cycle))
# ----------------------------------------------------------------------
# sample output
# 6->8->7->9->6->5->4->2->1->0->3->2->6
