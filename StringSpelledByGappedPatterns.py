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
    print ('deBruijn graph of gappedPatterns:', deBruijnGP)

# Find eulerian cycles and paths for debruijn paired graphs
#    eulerianPathFP = construct_eulerian_path(deBruijnFP)
#    eulerianPathSP = construct_eulerian_path(deBruijnSP)
#    print ('Eulerian path for First Patterns:', eulerianPathFP, '\nEulerian path for Second Patterns:', eulerianPathSP)

# Reconstruct strings from eulerian paths of first and second patterns
#    prefixString = reconstruct_string(eulerianPathFP)
#    suffixString = reconstruct_string(eulerianPathSP)
#    print ('Reconstructed string for first patterns:', prefixString, '\nReconstructed path for second patterns:', suffixString)



#    prefixEnder = kVal + dVal + 1
#    for i in range(prefixEnder, len(prefixString) + 1):
#        suffixOverlap = i - kVal - dVal
#        if prefixString[i] != suffixString[suffixOverlap]:
#            return 'There is no string spelled by the gapped patterns'
#        return prefixString + suffixString[-(prefixEnder - 1):]
    return "done"
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
#correctAnswer = 'GAAAGGTACAAATACTGGCGACCTCGCTGTTCGACACTTCATCACTGCTCCGGGGCGCTCAGGAGGGACGGTTCCCTGTACCATTGGAAGTCAATAGTCTAAGGTACAAAGAGAAGACCCGACCCGACAGAGGGGGTTCTGCGCCGGGTTTCGAGCTTGTAACCCCCCAGAGAATTAGATCCACCGTCTGTGTGGACAAAGTAGTAAAGCTAGCATACCAAATTGAAATTCGGAGTTTGACTACCAGATCCACGCATACGCTGCACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTAGAAATTCAGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGGGGTAATTCGTAGTTAGGTACAGAAAACTCCCGGACAGAACCGCATATAACCGATGAAGCAAGGGTTCTTCATTTAATACGACCCTAACCGGTATTGCTGCTAGCTTGATTTTCCTAGCAATCTAAACTCTATGTATGAGGCCACTCGGACGCCCGCTAGTGCCGGCAGCTAGCTACTGCCCTTCACCAGGAGCACGCACTATGCCTATCGGGCAATGCTGATCATACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATCCCTCTGCAGAAAGCGGTGGCGGCGGGTCTAAGCAAGTCCAACGCAATACCAGGAAATCACCGTATCGTTAGCGACCAGTAGGTGATGGTTTGTAAGTTCGGACTACAGGCGGATGTGTCCCCGCCAGTTAAAAGTCGACTTTCTGTTACAACTGCTCCCTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTCTAATGATCCCCACAGATCGTGTTTCAACGTTGAAGTCTGAATGGGTTCGTGAATATAGCCATCCAACGTGGACAATAAGATGAGCTTTATAGTTTCCGATCCTCATGGCGATCGAATAAGATCTATCCGCTTGTGTGTGTACGAGTCGCCGACTAACCGGTCTTGGGATATATACGTCACAGATTAAGTACTCGTCACGAGCTTGAATGGGAAGATAAGTAGACTCTTGTCGGGCACACACAGAGACTCCGACGCATCGAGATCGCAAACACTGCCTCCAGCCGGGGGATGCTAATCGTCGCGGTCGGTCCGAGCTTTATTCTACATCGTGGTGTTTCCGACCGAGCCATAAGAACAGTGTCCAAGTCACAAGAGGAGCACGCGGTGGAGGTCGTTCGCTATACAATATATTTGCAACTGTGTCTGGCATCACGCGCATTTCTCACACTTCCAAACGTGCTGCATTCTAAATGATTTCATGAATAGATTGTCTACTAGTTCACCCAAGGTATTACAGCACTGGTCATGTGCCGCTCTGGCACGGCTAGTATCAGGGCCGACTGTGTCCTAGGCGGCTGTTTTCGGGAGCCCAAGGGAAGCAATCAATGCGTTGACTGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATGACCCTAGGGCGGGCTGTTTGAGTGGGCTATCGGCGACCATTTCGCCCCGCAAGCCCCCGTCACGATACCGAGACCCTGAAGCTCATAAACGCCTATCTTTGTTGCATGAACAACGGGAGTAAGCGAGGCCAGGCCATACGGTTTCGGAGCCGCAGAATAGCCTTTACACGACCTCTACACAACCCAAAGTGAATATCCACGGGGTATTGTTTGTGCTGCTACCATGCGCAGTCAACATCGCCCACCGGCGATGTGTTTCAGATCTGCAGGCCCATCAACCGTTGTGACACCACCCCGGCTTTCAGAAGCAGTATGTCGGACACATTGACCTGTAGCGTTAGTTCTGTACAAGGGACCCTGCTCACTCGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTTGGGAATCTAATGCGGTCTGCCATGGGACCCCTAACAACAAGGGACGCCTCCACCTGTCTAGGAGGAAAGACTTTACACACATCTTCTAGTTTCGAGAAGCACCGTAGCCAGTGGACCCTGAGTGGTTACAAACAAATCGCAGTTTAGCGCTTACCGACAAAGGCGGGAGCTTCGTTCACATTAGGATTGAAAGAACTTAAGAGTCTGTAGGCTCGGAGGTCTCTATATATACCATCTAGTCGTCCGGGGCATGATTAACTAAGAGTTGATCTGAGTCGGAACATAATGCCTGATCTGACCCCAATTCACTACGGTCGCACGTTCCGGGAACACCTACCGATCAGTCCGGAACTGTGACCTAAGAAGTCTCAAGCCTTTACGTCAATGTTCCCGGTGAAGGACTGTGTAACGGTCGCCTTCGCGCCCCCCCATAGGCCGGTCCTTCTCGTTGCAGGATAGCTAAGTCCCATATAGAGTTGTTGGTGTACCATTACGCTGATTTTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGGTGTGTTCAACTCAGCATACCCGGTTAGTCTGGAGCACTCCCCGTGCCTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGCCGAAATCCTTGTCAAGTAGTGGCATCTACTGCCGCGGGGCAGGACTGATGCTGACCCAAGACCACGCTCCTATCAGCCGGTGGCGCATCAGGGTGGGCTATAAGTTATATTCCTACTGTACGGCTGAAGTCAGCTGTAGTCAGGGAGCGGTTCCTGAGCCGGCTGATTCCGCTCGTAATGCGCTATGTAGAAGCATAGTTAGCCTCGCGCTCGTGTGTGGGCAGTCGTATAAGTAGTTTAGCTCCCGATGCGAAGGAGTTGCAAGTACCTACAAACTCGTAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGTCCGGTTCATATTAACCATGCACCAAGGTTTGACTAAATCAACTCGTGGGAATCCGACGTGACAAAATCCCCAGATATGCCGGGGGTGCACGTGAATACGTCGTAAGTTGAGCGTCCTATGACGGGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGCCAAAGTAATACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGAATCCTTGACTGGCCTGCCGAGTGTTATCTCGCTGCTATCCCCCCCCAAGGTAGAAATGGAAGTGGGATCCAGCGCACCAAGGCACTTCACACAGGCATTACCCCAGCACCACGAATTAGCTTGCAGCTAAAGACAGGGTATTTTACGGAGTATATGATCTCTGTGAGGTACCGTATTCACACATCGTGGGATGTCTGCCATGAGCTTTTCCATTAGTATCCGGCGAGTTTTGATCCAAGTTACCAAACAAGGTTGTCCTCCAGGTCCTACGTGCTGAACGGCCAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGCCGGATCGGCCCTGACTCAAGTAAGGTCTGGTTGCTTGTCACTACATAAAGCCACGGAAGGGTGCGCGGCCCCAAAGCTGCGTCCGGATTCGACTCCCGTTTGCCTGGCTTCTCGACGAATCTAACGTTCTCATTAACCGAAAACCCTGAGCGGCTTAACCTCATTCGTCCCAGAATCAAACCCATCGTGTATCACCGTTGGCCCAGCAGGGAAGACAAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGTGGCTCCCCAACCGAAGAATAAGATCCCTTCGCCGCCACAGAAGCAACACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATTTAAAGTCTACCGTGGGGAGCCGGACGAGAACAAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGTAGGCGGCTACCTTCTGCCCATATCTCGGAATCCTTAGGGTCTTTAGCTGCTGGCAACACGGGGATGCCTTAGTCGCGGGGGAGCCGTTAGATCGGTTTCAGACTTGCCGACACCGTTCACCGTGGCCGGCCAACCGGCCGGTTCTGTGCTCCACTGAGTTTAGATAGGAATCCATACCTATAGTTTTACGTACACCAACTGGTTAACAAGCCGTCTCCGCGATGACAAACGGTGGGGGCACGAGCTCGAGGTAAGAGGTTGCGATCCGATTACAATGTATGACTACTTATAATGGTCTTACCTTAAATAGGGAACGGGTTTACAGATTACATGTCGGCGAAGACGTTACTGTATTTCGGCCAATGAGCAAATTCCCACCAGGCTCGCTCGGCTTAAATAGAATAGTCAGATTGGCTCTGAATGCTTTGGGCGTGCATTGAGAGCAGCCATATGTAATATTAAATGGCAGTACAAATCATACAGTTCAGAACTGCCGACAGCGCAGGAGTTTAAGGGTATCGAATATTGCGCTATCCGTGAGTGCTCTTAGCGATGGGGGGGCGGACCCTAAGTCTGACCCCCTCTCCTACCTTCTACGGATTACTATTATTGGCACTTGATGAGTAATCATTTCTAGCAAGAGTCTTATAAGGTAAACAATAACTTAGTAAGTGAGTCATGTAGTGTGCTTCCAGGACGAGTCGGCAAACTCTGTAGTCTTATGCTCATGTCTGACCTGCTGTGCCCAAAATTCTCTTCGTAAGGAGGGCTTTATAATGTTATGGGCACGACTTCGCATTGGGTCCACGCCCCAGGACTTCAGCATGTTATTTTGGGTTGCAGGATTTAAGAGAGCCTCATGCGTTGATAAGCCAAAGTGGGGGTATGGTGGGACCTCTCACCATGAGAGTTAAGTTAACTCACCGTGGCTCAAAAAAAGCTGGTTAGAATCTGCGAGTAATACGAGCGGGAAAATCTGGAATAACAGAAGCGACACCCTGACCTACAGTCGTTCAGTACTAGGTTACAAGTGAACCACTCGCGGATATAGTCAGGCGGGGATGTCCCGCGCGTTGATTAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGCGAGTTAATTGAAGTTTGCCTAGACACCGCTGAGGCTGGTTCGACATACCCTTAGGGAGGCCAAGCTATATAAAACCAAGATCATTGACCCCCTACGTGATACGTGATTTCAAACTTTACAATCATTAGGGTCGCCAGTGGAGAATCTATAGAATCTTTTCTACAGGCTACAGAGAAGCATTTTTCACAGGACCGCGTGGCGCAAACAATCCGATGGGGACCATCTGTGAACTCCCATACGTGACTATTCTGTGTCACATGAGGGGAGCTAGGGGGATTGAGTGCTCATGTCGGTTGGAGACCATTTTGAGTGCACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGCGAAGCATACTTACCTTGATCAACGCAGTGATTATTCATCTGAAGAGGATTGGGATAATACTGCGAACATATTGGAAAATTAACTGATTTATCTTCTGATCGATTCCCACACTCCACGAATTGGGGTGCCATGCTCCCATAGTAGGCCCTAGAGATGCCGATCATTCCGCAGGTGTGCCTAAGTGGACAGTCACTTGGCACTTAGGCCAATAAGTACAACAAAGGGATCAGTGGGCAAATTATCAGCGTACAATTCCCAGATATATAGGCGGCGAGAAAAGCTTCAAAAGACTTAATTTACTAGCCTCCTACAAACTCTAGATGAGGATTGGCTCTTGATGCTAGCGTTTTCATTTTCCATTACAAGACATTAGGCTGATAATTGCAGAGATTGGCGGCGTAGACTGACAGTCGCGATCAATCTGCGTGTTA'
#print ('Length of my answer:', len(myAnswer))
#print ('myanswer:', '\n', myAnswer)
#print ('\nLength of correct answer:', len(correctAnswer), '\n\ncorrectAnswer:', '\n', correctAnswer)
# -----testing------------------------------------------------------------
def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count
#print ('\n\nHammingDistance:', HammingDistance(myAnswer, correctAnswer))
