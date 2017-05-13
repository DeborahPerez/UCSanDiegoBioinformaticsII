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
#   VERSION:    20170512
########################################################################
# ****Find amino acid masses for each circularized petide option*****
# Summary
# *1. Experimental spectrum is "spectrum"
# *2. Collection of candidate linear peptides is "peptides"
# *3. Initialize "peptides" as empty
# *4. Expand "peptides" to contain all linear peptides of length 1
# *5. (Branching step) Continue step 4 process creating 18 new peptides
#    of length k + 1 for each amino acid string in peptide of length k
#    in "peptides" by appending all possible amino acid mass to the end
#    of peptide
# *6. (Bound step) For every expansion step of "peptides", trim "peptides"
#    by keeping only linear peptides consistent with experimental
#    "spectrum"
# *7. For every trim, check if new linear peptides have mass equal to
#    mass in the spectrum.
# 8. If mass of peptide equals mass of experimental parent, circularize
#    peptide and check if it is part of solution for cyclopeptide
#    sequencing problem.

# Pseudocode
# *CyclopeptideSequencing(Spectrum)
#   *Peptides ← a set containing only the empty peptide
#   *while Peptides is nonempty
#       *Peptides ← Expand(Peptides)
#       for each peptide Peptide in Peptides
#           if Mass(Peptide) = ParentMass(Spectrum)
#               if Cyclospectrum(Peptide) = Spectrum
#                   output Peptide
#               remove Peptide from Peptides
#           else if Peptide is not consistent with Spectrum
#               remove Peptide from Peptides
import sys
# ---bnb_cyc_pep_seq----------------------------------------------------
# Finds a peptide from a spectrum of masses using brute force
# @param spectrum of masses
# @return peptide
# ----------------------------------------------------------------------
def bnb_cyc_pep_seq(spectrum):
    aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                          'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                          'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                          'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
# Function to find a theoretical peptide match for experimental spectrum
    theoreticalMatch = _theorPepMatch(spectrum)
#    print (theoreticalMatch)

# Function to circularize chosen peptide and give spectrum values
    circularSequences = _circularize(theoreticalMatch)
    print (circularSequences)
    masterMasses = []
    for circularSequence in circularSequences:
        print (circularSequence)
        cirSeqMass = []
        for aminoAcid in circularSequence:
            mass = aminoAcidMassTable[aminoAcid]
            cirSeqMass.append(mass)
        masterMasses.append(cirSeqMass)
    return masterMasses


# ----------------------------------------------------------------------
# ---_theorPepMatch-----------------------------------------------------
# Finds a theoretical peptide match for experimental spectrum
# @param spectrum
# @return Theoretical peptide match
# ----------------------------------------------------------------------
def _theorPepMatch(spectrum):
    aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                          'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                          'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                          'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
    sortedSpectrum = sorted(spectrum)
    spectrumLength = len(spectrum)
    parentMass = sortedSpectrum[spectrumLength - 1]
    aminoAcidList = list(aminoAcidMassTable)
#    print ('Mass of peptide:', parentMass,
#    '\n\nsorted spectrum:', sortedSpectrum,
#    '\n\nlength of spectrum:', spectrumLength,
#    '\n\namino acids list:', aminoAcidList)
# Build empty peptides first
    peptides = []
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
# List of initialized candidate linear peptides
#    print ('\npeptides:', peptides)
    lenPeptides = len(peptides)
    while lenPeptides > 0:
# Bound step
        boundPeptides = []
        for pepInfo in peptides:
            pepList = pepInfo[0]
            stringPep = pepList[0]
            pepMassList = pepInfo[1]
            pepMass = pepMassList[0]
# Branch step
            if pepMass in spectrum:
                boundPeptides.append(pepInfo)
                if pepMass == parentMass:
                    if cyclic_spectrum(stringPep) == sortedSpectrum:
                        return stringPep
        peptides = expand_pep(boundPeptides)
# ---_circularize---------------------------------------------------------
# For a circular peptide of length k, finds every possible linear
# sequence of length k forward and backward
# @param Peptide
# @return All possible linear sequences
# ----------------------------------------------------------------------
def _circularize(peptide):
    lenPeptide = len(peptide)
    circularizedPeps = []
#    print ('\npeptide:', peptide)
    for i in range(lenPeptide):
        possiblePep = ''
        if i == 0:
            circularizedPeps.append(peptide)
            circularizedPeps.append(peptide[::-1])
        else:
            str1Part = peptide[:i]
            str2Part = peptide[i:]
            newCircPep = str2Part + str1Part
            circularizedPeps.append(newCircPep)
            circularizedPeps.append(newCircPep[::-1])
#            print ('\ni:', i,
#                '\npeptide string starting at ith index to end:',
#                peptide[i:],
#                '\nwrap around text:', peptide[:i],
#                '\nnew circularized peptide option:', newCircPep,
#                '\ncurrent circularized peptides:', circularizedPeps
#                )
    return circularizedPeps

# ----------------------------------------------------------------------
# ---expand_pep---------------------------------------------------------
# Expands possible peptide combination by 1 using all 18 amino acids
# @param Peptide list with their masses
# @return Expanded peptide list with their masses
# ----------------------------------------------------------------------
def expand_pep(peptideList):
    aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                          'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                          'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                          'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
    aminoAcidList = list(aminoAcidMassTable)
    newPeptideList = []
    for pepInfo in peptideList:
        pepList = pepInfo[0]
        stringPep = pepList[0]
        pepMassList = pepInfo[1]
        pepMass = pepMassList[0]

#        print ('\nPeptide information:', pepInfo,
#        '\nIndex 1 of Peptides:', pepList,
#        '\nPeptide as string:', stringPep,
#        '\nIndex 2 of Masses:', pepMassList,
#        '\nMass of Peptide:', pepMass)

        for aA in aminoAcidList:
            newPair = []
            newPep = []
            newStringPep = ''
            newStringPep += stringPep
            newStringPep += aA
            newPep.append(newStringPep)
            newMass = []
            newCalcMass = pepMass + aminoAcidMassTable[aA]
            newMass.append(newCalcMass)
            newPair.append(newPep)
            newPair.append(newMass)
            newPeptideList.append(newPair)

#            print ('new Peptide:', newStringPep,
#            '\nnew mass:', newCalcMass,
#            '\nnew pair info:', newPair)

    return newPeptideList

# ---cyclic_spectrum----------------------------------------------------
# Finds the cyclic spectrum of a peptide
# @param Peptide, Table of Amino Acids and their masses
# @return cyclic spectrum of the peptide
# ----------------------------------------------------------------------
def cyclic_spectrum(peptide):
    aminoAcidMassTable = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,
                          'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                          'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                          'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
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
possiblePeptide = bnb_cyc_pep_seq(spectrum)
print (possiblePeptide)
