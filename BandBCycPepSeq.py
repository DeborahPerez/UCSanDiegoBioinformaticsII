########################################################################
#   USAGE:
#       python3 BandBCycPepSeq.py
#   DESCRIPTION:
#       Grow candidate linear peptides whose theoretical spectra are
#       consistent with the experimental spectrum.
#   BIOINFORMATICS II GENOME SEQUENCING:
#       Input: Experimental mass spectrum
#       Output: Possible peptides
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20170427
########################################################################
import sys
# ---b_BCycPepSeq----------------------------------------------------
# Finds a peptide from a spectrum of masses using brute force
# @param spectrum of masses
# @return peptide
# ----------------------------------------------------------------------
aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                      'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                      'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                      'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
def bnb_cyc_pep_seq(spectrum, aminoAcidMassTable):
    sortedSpectrum = sorted(spectrum)
    spectrumLength = len(spectrum)
    parentMass = sortedSpectrum[spectrumLength - 1]
    aminoAcidList = list(aminoAcidMassTable)
    print ('Mass of peptide:', parentMass,
    '\n\nsorted spectrum:', sortedSpectrum,
    '\n\nlength of spectrum:', spectrumLength,
    '\n\namino acids list:', aminoAcidList)
# Build empty peptides first
    peptides = [[['empty'], [0]]]
# Add new amino acid and corresponding masses
    for aA in aminoAcidList:
        newPair = []
        newAA = []
        newMass = []
        newAA.append(aA)
        newMass.append(aminoAcidMassTable[aA])
        newPair.append(newAA)
        newPair.append(newMass)
        peptides.append(newPair)
    print (peptides)
# Print amino acid from peptides list
    for pepItem in peptides:
        pepItemAA = pepItem[0]
        print (pepItemAA)

# Find a way to modify  peptides list using expand_pep function

# ----------------------------------------------------------------------
# ---expand_pep---------------------------------------------------------
# Expands possible peptide combination by 1 using all 18 amino acids
# @param Peptide list with their masses
# @return Expanded peptide list with their masses
# ----------------------------------------------------------------------
def expand_pep(peptideList):
    aminoAcidList = list(aminoAcidMassTable)
    newPeptideList = []
    for peptide in peptideList:
        for aA in aminoAcidList:
            newPeptide = []
            newPeptide.append(peptide)
            newPeptide.append(aA)
    return newPeptideList
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
rawData = sys.stdin.read().splitlines()
rawSpectrum = rawData[0]
rawSpectrum = rawSpectrum.split()
spectrum = []
for item in rawSpectrum:
    number = int(item)
    spectrum.append(number)
possiblePeptide = bnb_cyc_pep_seq(spectrum, aminoAcidMassTable)
print (possiblePeptide)
