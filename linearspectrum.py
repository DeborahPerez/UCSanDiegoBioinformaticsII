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
    print (prefixMass)
    return aminoAcidMassTable
# ----------------------------------------------------------------------
# ---MAINCODE------------------------------------------------------------
peptide = 'NQEL'
linearSpectrum = linear_spectrum(peptide, aminoAcidMassTable)
print (linearSpectrum)
