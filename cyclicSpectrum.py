########################################################################
#   USAGE:
#       python3 cyclicSpectrum.py
#   DESCRIPTION:
#       Find the cyclic spectrum of a peptide.
#   BIOINFORMATICS II GENOME SEQUENCING:
#       Input: An amino acid string peptide.
#       Output: The cyclic spectrum of peptide.
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20170402
########################################################################
import sys
# ---cyclic_spectrum----------------------------------------------------
# Finds the cyclic spectrum of a peptide
# @param Peptide, Table of Amino Acids and their masses
# @return cyclic spectrum of the peptide
# ----------------------------------------------------------------------
aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                      'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                      'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                      'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
def cyclic_spectrum(peptide, aminoAcidMassTable):
    prefixMass = []
    prefixMass.append(0)
    aminoAcidList = list(aminoAcidMassTable)
    peptideLength = len(peptide)
#    print ('amino acid list:', aminoAcidList)
#    print ('\nprefix mass list:', prefixMass,
#     '\nlength of amino acid mass table:', len(aminoAcidMassTable),
#     '\npeptide sequence:', peptide,
#     '\n---------------------------------------------------------', '\n')
    for i in range(1, peptideLength + 1):
        currentPeptideIndex = i-1
        aminoAcid = peptide[currentPeptideIndex]
        aminoAcidMass = aminoAcidMassTable[aminoAcid]
        for j in range(len(aminoAcidList)):
            aminoAcidOfJ = aminoAcidList[j]
            aAMassOfJ = aminoAcidMassTable[aminoAcidOfJ]
            if aminoAcidOfJ == aminoAcid:
                prefixOfI = prefixMass[i-1] + aAMassOfJ
                prefixMass.append(prefixOfI)
#                print ('i:', i,
#                '\nCurrent peptide index:', currentPeptideIndex,
#                '\ni minus 1:', i-1,
#                '\n prefixMass of i - 1:', prefixMass[i-1],
#                '\nAmino acid:', aminoAcid, 'aminoAcidMass:', aminoAcidMass, '\nj:', j,
#                '\naminoAcidOfJ:', aminoAcidOfJ, '\nprefixOfI:', prefixOfI,
#                '\nCurrent prefixMass List:', prefixMass, '\n')
#    print ('---------------------------------------------------------', '\n')
    peptideMass = prefixMass[peptideLength]
#    print ('peptide mass:', peptideMass,
#         '\n---------------------------------------------------------', '\n')
    cyclicSpectrum = [0]
    for i in range(peptideLength):
        for j in range(i + 1, peptideLength + 1):
            jMinusI = prefixMass[j] - prefixMass[i]
            cyclicSpectrum.append(jMinusI)
#            print ('i:', i,
#            '\nj:', j,
#            '\nprefix mass of j minus prefix mass of i:', jMinusI,
#            '\n')
            if i > 0 and j < peptideLength:
                subPeptideMass = peptideMass - jMinusI
                cyclicSpectrum.append(subPeptideMass)
#                print ('mass of subpeptide j and i:', subPeptideMass)
    return sorted(cyclicSpectrum)
# ----------------------------------------------------------------------
# ---MAINCODE------------------------------------------------------------
rawData = sys.stdin.read().splitlines()
peptide = rawData[0]
cyclicSpectrum = cyclic_spectrum(peptide, aminoAcidMassTable)
print (" ".join(str(x) for x in cyclicSpectrum))
