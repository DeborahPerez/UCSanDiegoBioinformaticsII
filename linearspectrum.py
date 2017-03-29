########################################################################
#   USAGE:
#       python3 linearSpectrum.py
#   DESCRIPTION:
#       Find the linear spectrum of a peptide.
#   BIOINFORMATICS II GENOME SEQUENCING:
#       Input: An amino acid string peptide.
#       Output: The linear spectrum of peptide.
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20170324
########################################################################
import sys
# ---linear_spectrum----------------------------------------------------
# Finds the linear spectrum of a peptide
# @param Peptide, Table of Amino Acids and their masses
# @return Linear spectrum of the peptide
# ----------------------------------------------------------------------
aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                      'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                      'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                      'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
def linear_spectrum(peptide, aminoAcidMassTable):
    prefixMass = []
    prefixMass.append(0)
    aminoAcidList = list(aminoAcidMassTable)
    peptideLength = len(peptide)
    print ('amino acid list:', aminoAcidList)
    print ('\nprefix mass list:', prefixMass,
     '\nlength of amino acid mass table:', len(aminoAcidMassTable), '\n')
    for i in range(1, peptideLength + 1):
        currentPeptideIndex = i-1
        aminoAcid = peptide[currentPeptideIndex]
        aminoAcidMass = aminoAcidMassTable[aminoAcid]
        for j in range(len(aminoAcidList)):
            aminoAcidOfJ = aminoAcidList[j]
            if aminoAcidOfJ == aminoAcid:
                prefixOfI = prefixMass[i-1] + aminoAcidMass
                prefixMass.append(prefixOfI)
                print ('i:', i,
                '\nCurrent peptide index:', currentPeptideIndex,
                '\ni minus 1:', i-1,
                '\n prefixMass of i - 1:', prefixMass[i-1],
                '\nAmino acid:', aminoAcid, 'aminoAcidMass:', aminoAcidMass, '\nj:', j,
                '\naminoAcidOfJ:', aminoAcidOfJ, '\nprefixOfI:', prefixOfI,
                '\nCurrent prefixMass List:', prefixMass, '\n')
    linearSpectrum = [0]

#    print ('pseudocode left to complete:\npeptideMass ← PrefixMass(|Peptide|)',
#    '\nCyclicSpectrum ← a list consisting of the single integer 0',
#    '\n    for i ← 0 to |Peptide| − 1',
#    '\n        for j ← i + 1 to |Peptide|',
#    '\n            add PrefixMass(j) − PrefixMass(i) to CyclicSpectrum',
#    '\n            if i > 0 and j < |Peptide|',
#    '\n                add peptideMass - (PrefixMass(j) − PrefixMass(i)) to CyclicSpectrum',
#    '\n    return sorted list CyclicSpectrum')
    return prefixMass
# ----------------------------------------------------------------------
# ---MAINCODE------------------------------------------------------------
rawData = sys.stdin.read().splitlines()
peptide = rawData[0]
linearSpectrum = linear_spectrum(peptide, aminoAcidMassTable)
print (linearSpectrum)
