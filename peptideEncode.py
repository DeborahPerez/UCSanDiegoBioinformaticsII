########################################################################
#   USAGE:
#       python3 peptideEncode.py
#   DESCRIPTION:
#       Find substrings of a genome encoding a given amino acid
#       sequence
#   BIOINFORMATICS II GENOME SEQUENCING:
#       Input: A DNA string Text, an amino acid string Peptide, and
#       the array GeneticCode.
#       Output:  All substrings of Text encoding Peptide (if any such
#       substrings exist).
#-----------------------------------------------------------------------
#   CREATED BY: Deborah Perez
#   VERSION:    20170215
########################################################################
import sys
import urllib.request
# ---peptide_encode--------------------------------------------------------
# Constructs substrings of a genome encoding a given amino acid sequence
# @param dna string text, amino acid string text
# @return list of all substrings of Text encoding Peptide (if any such
# substrings exist).
# ----------------------------------------------------------------------

# Pseudocode:
# 1. reverse translate peptide to rna and find all possible rnas
# 2. put rnas into a list and duplicate dictionary
# 3a: loop through each letter of peptide - for each peptide
# 3b. start a string and record the first codon
# 3c. record next codon
# 4. search dna for rnas all 6 reading frames
# 5. for each rna found in dna add to new list of encoded peptides
# 6. include duplicates

def peptide_encode(dna, peptide):
    print ("dna:", dna, "\npeptide:", peptide, "\ncodonTable:", codonTable)

#reverse translate peptide letters to rna
    singlePeptideEncodes = []
    for i in range(len(peptide)):
        singleEncode = []
        aminoAcid = peptide[i]
        allRnas = iter(codonTable)
        for rnas in allRnas:
            if aminoAcid == codonTable[rnas]:
                singleEncode.append(rnas)
        singlePeptideEncodes.append(singleEncode)
#    print ('singlePeptideEncodes:', singlePeptideEncodes)

#run function to find all protein combinations as rna
    proteinCombinations = peptide_combos(singlePeptideEncodes)
#    print ('proteinCombinations Function Result:', proteinCombinations)
    for combination in proteinCombinations:
        print (combination)

#transcribe dna to rna and find all 6 reading frames
    #first 3 reading frames
    rf1Transcription = _transcription(dna)
    rf2Transcription = rf1Transcription[1:]
#    print (len(rf2Transcription)%3)
    cat = "CATSMEOW"
    print ('cat:', cat, '\nlength of cat with remainders:',
    (len(cat)%3), '\ncat without remainders:', cat[:-2],
    '\nlength of cat without remainders:', len(cat[:-2]))
#    rf2Transcription = rf2Transcription[:]
    rf3Transcription = rf1Transcription[2:]
    #next 3 reading frames from reverse complement read from 5' to 3'
    reverseComplement = reverse_complement(dna)
    rf4Transcription = _transcription(_reverse(reverseComplement))
    rf5Transcription = rf4Transcription[1:]
    rf6Transcription = rf3Transcription[2:]
    print ('\noriginal dna string:', dna, '\nlength of dna string:', len(dna),
    '\nreverse complement of dna:', reverse_complement(dna),
    '\nlength of reverse complement:', len(reverse_complement(dna)),
    '\nreading frame 1 transcription:', rf1Transcription,
    '\nlength of frame 1:', len(rf1Transcription),
    '\nreading frame 2 transcription:', rf2Transcription,
    '\nlength of frame 2:', len(rf2Transcription),
    '\nreading frame 3 transcription:', rf3Transcription,
    '\nlength of frame 3:', len(rf3Transcription),
    '\nreading frame 4 transcription:', rf4Transcription,
    '\nlength of frame 4:', len(rf4Transcription),
    '\nreading frame 5 transcription:', rf5Transcription,
    '\nlength of frame 5:', len(rf5Transcription),
    '\nreading frame 6 transcription:', rf6Transcription)
    '\nlength of frame 6:', len(rf6Transcription),


#(Framework) Begin list of peptide combinations with last peptide but reverseSPE


#Recursive Function
# ---peptide_combos-----------------------------------------------------
# Constructs all protein Combinations given a list of single rna peptide
# encode options
# @param reverse list of single rna peptide encodes
# @return list of combinations of protein
# ----------------------------------------------------------------------
def peptide_combos(singlePeptideEncodes):
#reverse the single peptides encodes list
    reverseSPE = singlePeptideEncodes[::-1]
    print ('\nreverseSPE List:', reverseSPE)
    proteinCombinations = []

#Initialize proteinCombinations with first index of reversed list
    currentSPE = reverseSPE[0]
    for SPE in currentSPE:
        proteinCombinations.append(SPE)
    del reverseSPE[0]
    print ('proteinCombinations from 1st index:', proteinCombinations)

#loop through reversed list and create combinations
    for rSPE in reverseSPE:
        currentSPE = rSPE
#        print ('\ncurrentSPE Group:', currentSPE, '\ncurrent proteinCombinations:', proteinCombinations)
        newPCombinations = []
        for proteinCombo in proteinCombinations:
            for i in range(len(currentSPE)):
                pComboString = ''
                pComboString += proteinCombo
#                print ('pComboString before concatenating with cSPEString:', pComboString)
                cSPEString = ''
                cSPEString += currentSPE[i]
                pComboString = cSPEString + pComboString
                newPCombinations.append(pComboString)
#                print ('i:', i, '\ncurrentSPE[i]:', currentSPE[i], '\ncSPEString:', cSPEString, '\npComboString after concatenation:', pComboString, '\nnewPCombinations after pComboString concatenation:', newPCombinations, '\ncurrent proteinCombinations:', proteinCombinations, '\n')
            proteinCombinations = list(newPCombinations)
#    print ('\nlength of protein combinations:', len(proteinCombinations), '\nprotein combinations:', proteinCombinations)
#        else:
    print ('length of protein combinations with duplicates:', len(proteinCombinations))
    list(set(proteinCombinations))
    print ('length of unique protein combinations:', len(proteinCombinations), '\nlength of one protein:', len(proteinCombinations[0]))
    return proteinCombinations







#Continue by distributing consecutive singlePeptides to framework, appending to the
#lists
#FIX HOW singlePeptideEncodes ARE APPENDING by reversing i and j order

#loop through reverseSPE list of encodes minus first instance


# ---translation--------------------------------------------------------
# Constructs an amino acid sequence by translating rna into amino acids
# @param rna string
# @return amino acid string
# ----------------------------------------------------------------------
def _translation(rna):
    codonList = [rna[i:i+3] for i in range(0, len(rna), 3)]
    print

    protein = ''
    for codon in codonList:
        if codon in codonTable:
            if codonTable[codon] == 'X':
                return protein
            protein += codonTable[codon]
    return protein
# ----------------------------------------------------------------------
# ---_transcription-----------------------------------------------------
# Transcribes dna into rna by replacing all 'T's with 'U's
# @param dna one or more lines of dna text
# @return a string as the rna transcript
# ----------------------------------------------------------------------
def _transcription(dna):
# Initialize rna output string
    transcript = ''
# Transcribe all 'T's found to be 'U's
    for line in dna:
        for nucleotide in line:
            if nucleotide != 'T':
                transcript += nucleotide
            else:
                transcript += 'U'
# Return string output
    return transcript
# ----------------------------------------------------------------------
# ---reverse_complement-------------------------------------------------
# Converts a strand of dna to its reverse complement
# @param dna one or more lines of dna text
# @return string of nucleotides as reverse complement
# ----------------------------------------------------------------------
def reverse_complement(dna):
# Initializes string to be returned for reverse complement
    dnaComplement = ''
# Reverses the order of the input string
    for line in dna:
        reverseDna = _reverse(line)
# Applies complements of reversed input to dnaComplement
        for nucleotide in reverseDna:
            dnaComplement += _complement(nucleotide)
# Returns output
    return dnaComplement
# ----------------------------------------------------------------------
# ---_reverse-----------------------------------------------------------
# Reverses text
# @param string of text
# @return string text in reverse order
# ----------------------------------------------------------------------
def _reverse(text):
# Initializes string to be returned for reverse order
    reversedText = ''
# Appends the index of text correlated to the value of count
    count = len(text) - 1
    for character in text:
        reversedText += text[count]
        count -=1
# Returns output
    return reversedText
# ----------------------------------------------------------------------
# ---_complement--------------------------------------------------------
# Translates complement of nucleotides: A=T, T=A, C=G, and G=C
# @param string of text
# @return string text with complement translation
# ----------------------------------------------------------------------
def _complement(nucleotide):
# Initializes string to be returned for translated complements
    complement = ''
# Translates a dna nucleotide to its dna complement nucleotide
    if nucleotide is 'A':
        complement += 'T'
    elif nucleotide is 'T':
        complement += 'A'
    elif nucleotide is 'C':
        complement += 'G'
    elif nucleotide is 'G':
        complement += 'C'
# Returns output
    return complement
#-----------------------------------------------------------------------
# ---MAINCODE-----------------------------------------------------------
# Open url and read byte information
byteText = urllib.request.urlopen(
"https://stepik.org/media/attachments/lessons/96/RNA_codon_table_1.txt"
).read()
# Convert bytes to string
strText = (byteText.decode("utf-8"))
# Save text line by line as a list
listText = strText.splitlines()

# Mark stop codons with an "X"
for i in range(len(listText)):
    if len(listText[i]) != 5:
        listText[i] += 'X'

# Create Codon Table
codonTable = {}
# Fill Codon Table with parsed information
for i in range(len(listText)):
    delimiter = ' '
    codon = listText[i].split(delimiter)
    for j in range(len(codon)):
        delimiter = ' '
        aminoAcid = codon[j].split(delimiter)
        codonTable[codon[0]] = codon[1]

# Standard Input
rawData = sys.stdin.read().splitlines()
dna = rawData[0]
peptide = rawData[1]
print (peptide_encode(dna, peptide))
