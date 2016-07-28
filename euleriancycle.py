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
#   VERSION:    20160727
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
    print (gMap)
# while loop for gMap with existing contents
    while len(gMap) > 0:
        print ("Gmap Length: ", len(gMap), "\nNode: ", node, "\nCycle: ", cycle, "\ngMap", gMap)
        if node in gMap:
            step = gMap[node]
            print ("\nStep: ", step)
#chooses first edge as next node aka step has more than 1 edge
            if len(step) > 1:
                cycle.append(step[0])
                edge = step[0]
                del gMap[node][0]
                node = (", ".join(edge))
                print ("\nNode: ", node)
            else:
                cycle.append(", ".join(step))
                edge = step
                del gMap[node]
                node = (", ".join(edge))
#In the case that gMap is not empty but node is not in gMap, we will traverse
#the cycle and choose an index existing in gMap and continue from there
        else:
            for i in range(len(cycle) - 1):
                print ("i:", i, "\nlength of cycle:", len(cycle), "\ncycle[i]", cycle[i])
                startPoint = cycle[i]
                if startPoint in gMap:
                    print ("i:", i, "\nstartPoint:", startPoint)
                    cutOffPoint = i
                    print ("cutOffPoint", cutOffPoint)
                    node = startPoint
#Before deleting part of cycle, find way to add keys and values back to dictionary
            del cycle[-int(cutOffPoint):]




    print (len(gMap), gMap)
    print ('answer', '6, 8, 7, 9, 6, 5, 4, 2, 1, 0, 3, 2, 6')
    return cycle

#
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
