########################################################################
#   USAGE:
#       python3 BFCyclopeptideSequencing.py
#   DESCRIPTION:
#       Find a possible peptide from an experimental mass spectrum
#       using a brute force algorithm
#   BIOINFORMATICS II GENOME SEQUENCING:
#       Input: Experimental mass spectrum
#       Output: peptide
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20170402
########################################################################
import sys
# ---BFCyclopep_seq----------------------------------------------------
# Finds a peptide from a spectrum of masses using brute force
# @param spectrum of masses
# @return peptide
# ----------------------------------------------------------------------
aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                      'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                      'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                      'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
def BFCyclopep_seq(spectrum, aminoAcidMassTable):
    sortedSpectrum = sorted(spectrum)
    spectrumLength = len(spectrum)
    parentMass = sortedSpectrum[spectrumLength - 1]
    for peptide in spectrum:
        if peptide == parentMass:
            if spectrum == cyclic_spectrum(peptide, aminoAcidMassTable):
                 return peptide
# ----------------------------------------------------------------------

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

    peptideMass = prefixMass[peptideLength]
    cyclicSpectrum = [0]
    for i in range(peptideLength):
        for j in range(i + 1, peptideLength + 1):
            jMinusI = prefixMass[j] - prefixMass[i]
            cyclicSpectrum.append(jMinusI)

            if i > 0 and j < peptideLength:
                subPeptideMass = peptideMass - jMinusI
                cyclicSpectrum.append(subPeptideMass)

    return sorted(cyclicSpectrum)
# ----------------------------------------------------------------------
# ---MAINCODE------------------------------------------------------------
#rawData = sys.stdin.read().splitlines()
#peptide = rawData[0]
spectrum = [0]
possiblePeptide = BFCyclopep_seq(spectrum, aminoAcidMassTable)
print (possiblePeptide)
