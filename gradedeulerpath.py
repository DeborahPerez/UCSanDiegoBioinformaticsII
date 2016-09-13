########################################################################
#   USAGE:
#       python3 Euleriancycle.py
#   DESCRIPTION:
#       Construct an Eulerian path from an adjacency list of an
#       Eulerian directed graph
#   BIOINFORMATICS II GENOME SEQUENCING:
#       1.2 From Euler's Theorem to an Algorithm for Finding Eulerian
#       Cycles - Step 6
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20160913
########################################################################
import sys
# ---construct_eulerian_path-------------------------------------------
# Constructs an Eulerian path from an adjacency list of an Eulerian
# directed graph
# @param Eulerian directed graph
# @return Eulerian path
# ----------------------------------------------------------------------
def construct_eulerian_path(gMap):
# Find Chooses a beginning node by finding the node that has more outgoing edges
# than incoming edges. Als
    cycle = []
    gMapKeys = list(gMap.keys())
    for listKey in gMapKeys:
        incomingEdges = 0
        outgoingEdges = gMap[listKey]
        for outgoingEdge in outgoingEdges:
            if outgoingEdge not in gMap:
                unbalancedNode = outgoingEdge
                unbalancedInOutNode = unbalancedNode
                gMap[unbalancedNode] = []
        for dictionaryKey in gMap:
            edges = gMap[dictionaryKey]
            if len(edges) != 1:
                for i in range(len(edges)):
                    if listKey == edges[i]:
                        incomingEdges += 1
            else:
                if listKey == edges[0]:
                    incomingEdges += 1
# Find unbalanced nodes
        if len(outgoingEdges) - incomingEdges == 1:
            unbalancedOutInNode = listKey
            activeNode = listKey
            cycle.append(activeNode)
        if incomingEdges - len(outgoingEdges) == 1:
            unbalancedInOutNode = listKey
            unbalancedNode = listKey

# Connect nodes in gMap
    gMap[unbalancedNode].append(activeNode)

# execute eulerian cycle
    eulerianCycle = construct_eulerian_cycle(gMap)

# Conversion of eulerian cycle to eulerian path requires rearrangement of
# euleriancycle

    for i in range(len(eulerianCycle)-1):
        if i == len(eulerianCycle):
            break
        if eulerianCycle[i] == unbalancedInOutNode and eulerianCycle[i+1] == unbalancedOutInNode:
            eulerianPath = eulerianCycle[i+1:] + eulerianCycle[1:i+1]
            return eulerianPath
# ----------------------------------------------------------------------

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
        outNodes = nodes[j].split(delimiter)
        gMap[nodes[0]] = outNodes

path = construct_eulerian_path(gMap)
print ("->".join(path))
# ----------------------------------------------------------------------
