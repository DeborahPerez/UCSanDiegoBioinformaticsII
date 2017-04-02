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
# Finds a peptide from a spectrum of masses
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
    peptideMass = sortedSpectrum[spectrumLength - 1]
    return 'true'
# ----------------------------------------------------------------------
# ---MAINCODE------------------------------------------------------------
#rawData = sys.stdin.read().splitlines()
#peptide = rawData[0]
spectrum = [0]
possiblePeptide = BFCyclopep_seq(spectrum, aminoAcidMassTable)
print (possiblePeptide)
