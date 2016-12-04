########################################################################
#   USAGE:
#       python3 StringSpelledByGappedPatterns.py
#   DESCRIPTION:
#       Implement StringSpelledByGappedPatterns.
#   BIOINFORMATICS II GENOME SEQUENCING:
#       Input: Integers k and d followed by a sequence of (k, d)-mers
#       (a1|b1),… , (an|bn) such that Suffix(ai|bi) = Prefix(ai+1|bi+1)
#       for 1 ≤ i ≤ n-1.
#       Output: A string Text of length k + d + k + n - 1 such that the
#       i-th (k, d)-mer in Text is equal to (ai|bi)  for 1 ≤ i ≤ n (if
#       such a string exists).
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20161118
########################################################################
import sys
# ***1. STRING RECONSTRUCTION FOR SMALL DATASET
#    2. STRING RECONSTRUCTION FOR LARGE DATASET
#    3. Check large dataset answer against homework key


# ---string_spelled_by_g_p----------------------------------------------
# Reconstructs a string by overlapping a prefixString and suffixString
# through firstPatterns and secondPatterns from gapped patterns
# @param list of gapped patterns, k, and d values
# @return reconstructed contiguous string
# ----------------------------------------------------------------------
def string_spelled_by_g_p(gappedPatterns, kVal, dVal):
    firstPatterns = []
    secondPatterns = []
    for kDmers in gappedPatterns:
        firstPatterns.append(kDmers[0])
        secondPatterns.append(kDmers[1])
#    print ('firstPatterns:', firstPatterns, '\n\nsecondPatterns:', secondPatterns)

#    prefixString = reconstruct_string(firstPatterns)
#    suffixString = reconstruct_string(secondPatterns)
#    print ('prefixString:', prefixString, '\n', 'suffixString:', suffixString, '\n')

# Find deBruijn graph for gappedPatterns
    deBruijnGP = paired_kmer_debruijn(gappedPatterns)
#    deBruijnSP = kmer_debruijn_graph(secondPatterns)
##    print ('deBruijn graph of gappedPatterns:', deBruijnGP)

# Find eulerian cycles and paths for debruijn paired kmers graph
#    eulerianPathFP = construct_eulerian_path(deBruijnFP)
#    eulerianPathSP = construct_eulerian_path(deBruijnSP)
#    print ('Eulerian path for First Patterns:', eulerianPathFP, '\nEulerian path for Second Patterns:', eulerianPathSP)
    eulerianPathRP = construct_eulerian_path(deBruijnGP)
##    print (eulerianPathRP)
# Reconstruct strings from eulerian paths of first and second patterns
#    prefixString = reconstruct_string(eulerianPathFP)
#    suffixString = reconstruct_string(eulerianPathSP)
#    print ('Reconstructed string for first patterns:', prefixString, '\nReconstructed path for second patterns:', suffixString)
    fPatterns = []
    sPatterns = []
    for step in eulerianPathRP:
        fPatterns.append(step.split('|')[0])
        sPatterns.append(step.split('|')[1])
#    print ('firstPatterns:', fPatterns, '\n\nsecondpatterns:', sPatterns)

    prefixString = reconstruct_string(fPatterns)
    suffixString = reconstruct_string(sPatterns)
##    print ('\nprefixString:', prefixString, '\nsuffixString:', suffixString, '\n')

    prefixEnder = kVal + dVal + 1
##    print ('prefixEnder:', prefixEnder, '\nlength of prefix string:', len(prefixString) + 1)
    for i in range(prefixEnder, len(prefixString) + 1):
        suffixOverlap = i - kVal - dVal
##        print ('suffix overlap:', suffixOverlap, '\ni in prefixString:', prefixString[i], '\nsuffix string overlap in the suffix string:', suffixString[suffixOverlap])
        if prefixString[i] != suffixString[suffixOverlap]:
            return 'There is no string spelled by the gapped patterns'
        return prefixString + suffixString[-(prefixEnder - 1):]
#    return "done"
# ----------------------------------------------------------------------

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
##    print ('\nall gmap keys:', gMapKeys)
    for listKey in gMapKeys:
        incomingEdges = 0
        outgoingEdges = gMap[listKey]
##        print ('\n\nlistKey(key in list of keys we are looping):', listKey)
##        print ('----ListValue', outgoingEdges)
#        for outgoingEdge in outgoingEdges:
        if outgoingEdges not in gMap:
##            print ('outgoing edge not in gMap!:', outgoingEdges)
            unbalancedNode = outgoingEdges
            unbalancedInOutNode = unbalancedNode
            gMap[unbalancedNode] = []
##            print ('****gMap after addition of unbalanced edge:', '\n', gMap)
        for dictionaryKey in gMap:
            edges = gMap[dictionaryKey]
#            print ('\ndictionaryKey:', dictionaryKey, '\nvalues:', edges)
            if listKey == edges:
                incomingEdges += 1
##        print ('incomingEdges:', incomingEdges, '\noutgoingEdges:', outgoingEdges, '\n-------------------------------------------------------')
        outgoingEdges = outgoingEdges.count('|')
        outMinusIn = outgoingEdges - incomingEdges
        inMinusOut = incomingEdges - outgoingEdges
##        print ('outgoing minus incoming:', outMinusIn)
##        print ('incoming minus outgoing:', inMinusOut)
#            if len(edges) != 1:
#                for i in range(len(edges)):
#                    if listKey == edges[i]:
#                        incomingEdges += 1
#            else:
#                if listKey == edges[0]:
#                    incomingEdges += 1
# Find unbalanced nodes
#        if len(outgoingEdges) - incomingEdges == 1:
        if outMinusIn == 1:
##            print ('new active node')
            unbalancedOutInNode = listKey
            activeNode = listKey
            cycle.append(activeNode)
        elif inMinusOut == 1:
##            print ('MEMEMEMEMEM')
            unbalancedInOutNode = listKey
            unbalancedNode = listKey

# Connect nodes in gMap
    gMap[unbalancedInOutNode] = unbalancedOutInNode

# execute eulerian cycle
    eulerianCycle = construct_ec_from_pr(gMap)

# Conversion of eulerian cycle to eulerian path requires rearrangement of
# euleriancycle

    for i in range(len(eulerianCycle)-1):
        if i == len(eulerianCycle):
            break
        if eulerianCycle[i] == unbalancedInOutNode and eulerianCycle[i+1] == unbalancedOutInNode:
            eulerianPath = eulerianCycle[i+1:] + eulerianCycle[1:i+1]
            return eulerianPath
# ----------------------------------------------------------------------

# ---construct_ec_from_pr-----------------------------------------------
# Constructs an Eulerian cycle from an adjacency list of paired reads
# an Eulerian directed graph
# @param Eulerian directed graph of paired reads
# @return Eulerian cycle of paired reads
# ----------------------------------------------------------------------
def construct_ec_from_pr(gMap):
##    print ('\n\n\nEulerian cycle code begins here*********')
    node = next(iter(gMap))
##    print ('first node:', node)
    cycle = []
    cycle.append(node)
##    print ('beginning of cycle:', cycle)

# while loop for gMap with existing contents
    while len(gMap) > 0:
        if node in gMap:
            step = gMap[node]
            cycle.append(step)
#dictionary deletion mechanism
#            if len(steps) > 1:
#                del gMap[node][0]
#            else:
#                del gMap[node]
            del gMap[node]
##            print ('\ncycle:', cycle, '\ngMap:', gMap)
            node = step
#            print ('cycle:', cycle)
#In the case that gMap is not empty and node is not in gMap, we will traverse
#the cycle and choose an index existing in gMap and continue from there
        else:
##            print ('node not in gMap:', node, '\ncurrent cycle:', cycle, '\ncurrent gMap', gMap)
            for i in range(len(cycle) - 1):
                startPoint = cycle[i]
##                print ('i:', i, 'startPoint:', startPoint)
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

# ---kmer_debruijn_graph------------------------------------------------
# Constructs a de Bruijn graph from a set of paired k-mers
# @param list of paired kmers
# @return deBruijn graph of paired kmers
# ----------------------------------------------------------------------
def paired_kmer_debruijn(pairedKmers):
    graphMap = {}
    prefixPRList = []
    suffixPRList = []
    for pairedKmer in pairedKmers:
        prefixPR = ''
        suffixPR = ''
#        print ('pairedKmer:', pairedKmer)
        for i in range (len(pairedKmer)):
            prefixPR += (_prefix(pairedKmer[i]))
            suffixPR += (_suffix(pairedKmer[i]))
            if i == 0:
                prefixPR += '|'
                suffixPR += '|'
        graphMap[prefixPR] = suffixPR
#    graphMap[prefixPR].append(suffixPR)
#        print ('prefixes:', prefixPR, '\nsuffixes:', suffixPR)
#        print ('\n')
#        prefixPRList.append(prefixPR)
#        suffixPRList.append(suffixPR)

#    print ("all prefix paired reads:", prefixPRList, '\n\n\n')
#    print ("all suffix paired reads:", suffixPRList)

#    for prefixPR in prefixPRList:
#        if prefixPR not in graphMap:
#            graphMap[prefixPR] = [_suffix(kmer)]
#        else:
#            graphMap[_prefix(kmer)].append(_suffix(kmer))
#    return "done"
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
rawData = sys.stdin.read().splitlines()
kVal, dVal = map(int, rawData[0].split(" "))
#print ('kVal=', kVal, '\ndVal=', dVal)
del rawData[0]

gappedPatterns = []
for line in rawData:
    gappedPattern = line.split("|")
#    print (gappedPattern)
    gappedPatterns.append(gappedPattern)
#print ('\ngappedPatterns:', '\n', gappedPatterns, '\n')

myAnswer = string_spelled_by_g_p(gappedPatterns, kVal, dVal)
#correctAnswer = 'GTGGTCGTGAGATGTTGA'
correctAnswer = 'GAAAGGTACAAATACTGGCGACCTCGCTGTTCGACACTTCATCACTGCTCCGGGGCGCTCAGGAGGGACGGTTCCCTGTACCATTGGAAGTCAATAGTCTAAGGTACAAAGAGAAGACCCGACCCGACAGAGGGGGTTCTGCGCCGGGTTTCGAGCTTGTAACCCCCCAGAGAATTAGATCCACCGTCTGTGTGGACAAAGTAGTAAAGCTAGCATACCAAATTGAAATTCGGAGTTTGACTACCAGATCCACGCATACGCTGCACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTAGAAATTCAGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGGGGTAATTCGTAGTTAGGTACAGAAAACTCCCGGACAGAACCGCATATAACCGATGAAGCAAGGGTTCTTCATTTAATACGACCCTAACCGGTATTGCTGCTAGCTTGATTTTCCTAGCAATCTAAACTCTATGTATGAGGCCACTCGGACGCCCGCTAGTGCCGGCAGCTAGCTACTGCCCTTCACCAGGAGCACGCACTATGCCTATCGGGCAATGCTGATCATACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATCCCTCTGCAGAAAGCGGTGGCGGCGGGTCTAAGCAAGTCCAACGCAATACCAGGAAATCACCGTATCGTTAGCGACCAGTAGGTGATGGTTTGTAAGTTCGGACTACAGGCGGATGTGTCCCCGCCAGTTAAAAGTCGACTTTCTGTTACAACTGCTCCCTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTCTAATGATCCCCACAGATCGTGTTTCAACGTTGAAGTCTGAATGGGTTCGTGAATATAGCCATCCAACGTGGACAATAAGATGAGCTTTATAGTTTCCGATCCTCATGGCGATCGAATAAGATCTATCCGCTTGTGTGTGTACGAGTCGCCGACTAACCGGTCTTGGGATATATACGTCACAGATTAAGTACTCGTCACGAGCTTGAATGGGAAGATAAGTAGACTCTTGTCGGGCACACACAGAGACTCCGACGCATCGAGATCGCAAACACTGCCTCCAGCCGGGGGATGCTAATCGTCGCGGTCGGTCCGAGCTTTATTCTACATCGTGGTGTTTCCGACCGAGCCATAAGAACAGTGTCCAAGTCACAAGAGGAGCACGCGGTGGAGGTCGTTCGCTATACAATATATTTGCAACTGTGTCTGGCATCACGCGCATTTCTCACACTTCCAAACGTGCTGCATTCTAAATGATTTCATGAATAGATTGTCTACTAGTTCACCCAAGGTATTACAGCACTGGTCATGTGCCGCTCTGGCACGGCTAGTATCAGGGCCGACTGTGTCCTAGGCGGCTGTTTTCGGGAGCCCAAGGGAAGCAATCAATGCGTTGACTGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATGACCCTAGGGCGGGCTGTTTGAGTGGGCTATCGGCGACCATTTCGCCCCGCAAGCCCCCGTCACGATACCGAGACCCTGAAGCTCATAAACGCCTATCTTTGTTGCATGAACAACGGGAGTAAGCGAGGCCAGGCCATACGGTTTCGGAGCCGCAGAATAGCCTTTACACGACCTCTACACAACCCAAAGTGAATATCCACGGGGTATTGTTTGTGCTGCTACCATGCGCAGTCAACATCGCCCACCGGCGATGTGTTTCAGATCTGCAGGCCCATCAACCGTTGTGACACCACCCCGGCTTTCAGAAGCAGTATGTCGGACACATTGACCTGTAGCGTTAGTTCTGTACAAGGGACCCTGCTCACTCGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTTGGGAATCTAATGCGGTCTGCCATGGGACCCCTAACAACAAGGGACGCCTCCACCTGTCTAGGAGGAAAGACTTTACACACATCTTCTAGTTTCGAGAAGCACCGTAGCCAGTGGACCCTGAGTGGTTACAAACAAATCGCAGTTTAGCGCTTACCGACAAAGGCGGGAGCTTCGTTCACATTAGGATTGAAAGAACTTAAGAGTCTGTAGGCTCGGAGGTCTCTATATATACCATCTAGTCGTCCGGGGCATGATTAACTAAGAGTTGATCTGAGTCGGAACATAATGCCTGATCTGACCCCAATTCACTACGGTCGCACGTTCCGGGAACACCTACCGATCAGTCCGGAACTGTGACCTAAGAAGTCTCAAGCCTTTACGTCAATGTTCCCGGTGAAGGACTGTGTAACGGTCGCCTTCGCGCCCCCCCATAGGCCGGTCCTTCTCGTTGCAGGATAGCTAAGTCCCATATAGAGTTGTTGGTGTACCATTACGCTGATTTTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGGTGTGTTCAACTCAGCATACCCGGTTAGTCTGGAGCACTCCCCGTGCCTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGCCGAAATCCTTGTCAAGTAGTGGCATCTACTGCCGCGGGGCAGGACTGATGCTGACCCAAGACCACGCTCCTATCAGCCGGTGGCGCATCAGGGTGGGCTATAAGTTATATTCCTACTGTACGGCTGAAGTCAGCTGTAGTCAGGGAGCGGTTCCTGAGCCGGCTGATTCCGCTCGTAATGCGCTATGTAGAAGCATAGTTAGCCTCGCGCTCGTGTGTGGGCAGTCGTATAAGTAGTTTAGCTCCCGATGCGAAGGAGTTGCAAGTACCTACAAACTCGTAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGTCCGGTTCATATTAACCATGCACCAAGGTTTGACTAAATCAACTCGTGGGAATCCGACGTGACAAAATCCCCAGATATGCCGGGGGTGCACGTGAATACGTCGTAAGTTGAGCGTCCTATGACGGGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGCCAAAGTAATACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGAATCCTTGACTGGCCTGCCGAGTGTTATCTCGCTGCTATCCCCCCCCAAGGTAGAAATGGAAGTGGGATCCAGCGCACCAAGGCACTTCACACAGGCATTACCCCAGCACCACGAATTAGCTTGCAGCTAAAGACAGGGTATTTTACGGAGTATATGATCTCTGTGAGGTACCGTATTCACACATCGTGGGATGTCTGCCATGAGCTTTTCCATTAGTATCCGGCGAGTTTTGATCCAAGTTACCAAACAAGGTTGTCCTCCAGGTCCTACGTGCTGAACGGCCAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGCCGGATCGGCCCTGACTCAAGTAAGGTCTGGTTGCTTGTCACTACATAAAGCCACGGAAGGGTGCGCGGCCCCAAAGCTGCGTCCGGATTCGACTCCCGTTTGCCTGGCTTCTCGACGAATCTAACGTTCTCATTAACCGAAAACCCTGAGCGGCTTAACCTCATTCGTCCCAGAATCAAACCCATCGTGTATCACCGTTGGCCCAGCAGGGAAGACAAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGTGGCTCCCCAACCGAAGAATAAGATCCCTTCGCCGCCACAGAAGCAACACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATTTAAAGTCTACCGTGGGGAGCCGGACGAGAACAAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGTAGGCGGCTACCTTCTGCCCATATCTCGGAATCCTTAGGGTCTTTAGCTGCTGGCAACACGGGGATGCCTTAGTCGCGGGGGAGCCGTTAGATCGGTTTCAGACTTGCCGACACCGTTCACCGTGGCCGGCCAACCGGCCGGTTCTGTGCTCCACTGAGTTTAGATAGGAATCCATACCTATAGTTTTACGTACACCAACTGGTTAACAAGCCGTCTCCGCGATGACAAACGGTGGGGGCACGAGCTCGAGGTAAGAGGTTGCGATCCGATTACAATGTATGACTACTTATAATGGTCTTACCTTAAATAGGGAACGGGTTTACAGATTACATGTCGGCGAAGACGTTACTGTATTTCGGCCAATGAGCAAATTCCCACCAGGCTCGCTCGGCTTAAATAGAATAGTCAGATTGGCTCTGAATGCTTTGGGCGTGCATTGAGAGCAGCCATATGTAATATTAAATGGCAGTACAAATCATACAGTTCAGAACTGCCGACAGCGCAGGAGTTTAAGGGTATCGAATATTGCGCTATCCGTGAGTGCTCTTAGCGATGGGGGGGCGGACCCTAAGTCTGACCCCCTCTCCTACCTTCTACGGATTACTATTATTGGCACTTGATGAGTAATCATTTCTAGCAAGAGTCTTATAAGGTAAACAATAACTTAGTAAGTGAGTCATGTAGTGTGCTTCCAGGACGAGTCGGCAAACTCTGTAGTCTTATGCTCATGTCTGACCTGCTGTGCCCAAAATTCTCTTCGTAAGGAGGGCTTTATAATGTTATGGGCACGACTTCGCATTGGGTCCACGCCCCAGGACTTCAGCATGTTATTTTGGGTTGCAGGATTTAAGAGAGCCTCATGCGTTGATAAGCCAAAGTGGGGGTATGGTGGGACCTCTCACCATGAGAGTTAAGTTAACTCACCGTGGCTCAAAAAAAGCTGGTTAGAATCTGCGAGTAATACGAGCGGGAAAATCTGGAATAACAGAAGCGACACCCTGACCTACAGTCGTTCAGTACTAGGTTACAAGTGAACCACTCGCGGATATAGTCAGGCGGGGATGTCCCGCGCGTTGATTAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGCGAGTTAATTGAAGTTTGCCTAGACACCGCTGAGGCTGGTTCGACATACCCTTAGGGAGGCCAAGCTATATAAAACCAAGATCATTGACCCCCTACGTGATACGTGATTTCAAACTTTACAATCATTAGGGTCGCCAGTGGAGAATCTATAGAATCTTTTCTACAGGCTACAGAGAAGCATTTTTCACAGGACCGCGTGGCGCAAACAATCCGATGGGGACCATCTGTGAACTCCCATACGTGACTATTCTGTGTCACATGAGGGGAGCTAGGGGGATTGAGTGCTCATGTCGGTTGGAGACCATTTTGAGTGCACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGCGAAGCATACTTACCTTGATCAACGCAGTGATTATTCATCTGAAGAGGATTGGGATAATACTGCGAACATATTGGAAAATTAACTGATTTATCTTCTGATCGATTCCCACACTCCACGAATTGGGGTGCCATGCTCCCATAGTAGGCCCTAGAGATGCCGATCATTCCGCAGGTGTGCCTAAGTGGACAGTCACTTGGCACTTAGGCCAATAAGTACAACAAAGGGATCAGTGGGCAAATTATCAGCGTACAATTCCCAGATATATAGGCGGCGAGAAAAGCTTCAAAAGACTTAATTTACTAGCCTCCTACAAACTCTAGATGAGGATTGGCTCTTGATGCTAGCGTTTTCATTTTCCATTACAAGACATTAGGCTGATAATTGCAGAGATTGGCGGCGTAGACTGACAGTCGCGATCAATCTGCGTGTTA'
#print ('Length of my answer:', len(myAnswer))
print ('myanswer:', '\n', myAnswer)
#print ('\nLength of correct answer:', len(correctAnswer), '\ncorrectAnswer:', '\n', correctAnswer)
# -----testing------------------------------------------------------------
def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count
#print ('\n\nHammingDistance:', HammingDistance(myAnswer, correctAnswer))
