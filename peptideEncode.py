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
#   VERSION:    20170323
#   LINE If location is len peptideencode - peptide length, then......
########################################################################
import sys
#import urllib.request
# ---peptide_encode--------------------------------------------------------
# Constructs substrings of a genome encoding a given amino acid sequence
# @param dna string text, amino acid string text
# @return list of all substrings of Text encoding Peptide (if any such
# substrings exist).
# ----------------------------------------------------------------------
def peptide_encode(dna, peptide):

    print ("\ndna:", dna, '\ndata type:', type(dna), "\npeptide:", peptide,
     "\n\ncodonTable:", codonTable, '\n')

#run function to find all protein combinations as rna
#    proteinCombinations = peptide_combos(peptide)
#    print ('proteinCombinations Function Result:', proteinCombinations)i
#    print ('\nall possible protein encode combinations:')
#    for combination in proteinCombinations:
#        print (combination)

#Create a list of lists(dna, translation, transcription) to keep track of
#translations and transcriptions
    conversionList = []

#find all 6 dna reading frames
    print ('\nAll 6 Reading frames:')
    sixDnaRFs = _findAllSixRFs(dna)
    for rf in sixDnaRFs:
        print (rf)
#    print ('\nall 6 dna reading frames:')
#    for dnaRF in sixDnaRFs:
#        print (dnaRF)

#loop to create a list of the transcriptions and translations
    for i in range(len(sixDnaRFs)):
        singleConversion = []
        singleConversion.append(sixDnaRFs[i])
#transcribe reading frame
#        print ('\ntranscribed reading frame:')
        transcribedRF = _transcription(sixDnaRFs[i])
#        print (transcribedRF, '\nlength of transcribedRF:', len(transcribedRF))
        singleConversion.append(transcribedRF)
#translate reading frame
        translatedRF = _translation(transcribedRF)
#        print ('\ntranslated frame:', translatedRF, '\nlength of frame:',
#        len(translatedRF))
        singleConversion.append(translatedRF)
        conversionList.append(singleConversion)
#conversion list now index zero is dna, index 1 is rna, and index 2 is protein
#    print (ConversionList)
#    for i in range(len(conversionList)):
#        print ('\nreading frame:', i + 1)
#        for j in range(len(conversionList[i])):
#            if j == 0:
#                print ('\ndna:', conversionList[i][j])
#            elif j == 1:
#                print ('\nrna:', conversionList[i][j])
#            elif j == 2:
#                print ('\nprotein sequence:', conversionList[i][j])

#Test on splicing to adjust reading frame length
#    cat = "CATSMEOW"
#    print ('cat:', cat, '\nlength of cat with remainders divisible by 3:',
#    (len(cat)%3), '\ncat divisible by 3:', cat[:3],
#    '\nlength of cat divisible by 3:', len(cat[:-2]))
#------------------------------------------------------------------------------

#Find peptide locations in first three reading frames by index 2
    subStringEncPept = []
    peptideLength = len(peptide)
    dnaPeptideLength = peptideLength * 3

    print ('\n----------------------------------------------------------',
    '\n*Find peptides in all reading frames*', '\npeptide:', peptide,
    '\nLength of peptide:', peptideLength,
    '\nLength of dna encoding peptide:', dnaPeptideLength,
    '\n----------------------------------------------------------', '\n')

    for i in range(3):
        dnaString = conversionList[i][0]
        peptideString = conversionList[i][2]
        rnaString = conversionList[i][1]
        proteinLocation = _locateSubtext(peptide, peptideString)

        print ('Value of i:', i, '\nReading frame:', i + 1,
        '\nReading frame', i+1, '\npeptide encode:', peptideString,
        '\nReading frame', i+1, 'dna string:', dnaString,
         '\nprotein location(s):', proteinLocation, '\n')

        if len(proteinLocation) > 0:
            for location in proteinLocation:
                ssPair = []
                startPointPep = location

                #starting point setting
                if location == 0:
                    startPointDna = 0
                else:
                    startPointDna = location * 3
                print ('\nStarting Point in DNA:',
                startPointDna,
                '\nStarting Point in DNA + dna peptide length:',
                startPointDna + dnaPeptideLength)

                #ending point setting
                endPointDna = startPointDna + dnaPeptideLength
                if endPointDna == len(dnaString):
                    subString = dnaString[startPointDna:endPointDna - 1]
                else:
                    subString = dnaString[startPointDna:endPointDna]

                ssPair.append(i+1)
                ssPair.append(subString)
                subStringEncPept.append(ssPair)

                print ('\nStarting point index in peptide and including:',
                startPointPep,
                '\nStarting point value in peptide and including:',
                peptideString[startPointPep],
#                '\nEnd of location index in peptide and including:',
#                endPointPep,
#                '\nEnd of location value in peptide and including:',
#                peptideString[endPointPep],
                '\nStarting point index in dna and including:', startPointDna,
                '\nStarting point value in dna and including:',
                dnaString[startPointDna],
                '\nEnd of location index in dna NOT including:', endPointDna,
#                '\nEnd of location value in dna NOT including:',
#                dnaString[endPointDna],
                '\nLength of dnaString:', len(dnaString),
                '\nLength of translation:', len(peptideString),
                '\n\nDna substring:', subString,
                '\nLength of substring:', len(subString), '\n')
#find peptide locations in reverse complement reading frames by index 2
    print ('\n')
    for i in range(3,6):
        dnaString = conversionList[i][0]
        peptideString = conversionList[i][2]
        rnaString = conversionList[i][1]
        proteinLocation = _locateSubtext(peptide, peptideString)

        print ('Value of i:', i,
        '\nReading frame:', i + 1,
#        '\nReading frame', i+1, 'peptide encode:', peptideString,
#        '\nReading frame', i+1, 'dna string:', dnaString,
         '\nprotein location(s):', proteinLocation, '\n')

        if len(proteinLocation) > 0:
            for location in proteinLocation:
                ssPair = []

                #starting point setting
                startPointPep = location
                if location == 0:
                    startPointDna = 0
                else:
                    startPointDna = location * 3
                print ('\nStarting Point in DNA:',
                startPointDna,
                '\nStarting Point in DNA + dna peptide length:',
                startPointDna + dnaPeptideLength,
                '\nlast character:', dnaString[startPointDna + dnaPeptideLength])

                #ending point setting
                endPointDna = startPointDna + dnaPeptideLength
                if endPointDna == len(dnaString):
                    subString = dnaString[startPointDna:endPointDna - 1]
                else:
                    subString = dnaString[startPointDna:endPointDna]

                subString = reverse_complement(subString)
                ssPair.append(i+1)
                ssPair.append(subString)
                subStringEncPept.append(ssPair)

                print ('\nStarting point index in peptide and including:',
                startPointPep,
                '\nStarting point value in peptide and including:',
                peptideString[startPointPep],
                '\nEnd of location index in peptide and including:',
                location,
                '\nEnd of location value in peptide and including:',
                peptideString[location],
                '\nStarting point index in dna and including:', startPointDna,
                '\nStarting point value in dna and including:',
                dnaString[startPointDna],
                '\nEnd of location index in dna NOT including:', endPointDna,
#                '\nEnd of location value in dna NOT including:',
#                dnaString[endPointDna],
                '\nLength of dnaString:', len(dnaString),
                '\nLength of translation:', len(peptideString),
                '\n\nDna substring:', subString,
                '\nLength of substring:', len(subString), '\n')

    return subStringEncPept
#    for subscript in subStringEncPept:
#        print (subscript)
# ---_locateDnaMotif----------------------------------------------------
# finds location of subtext within a text
# @param peptide/subtext and dnaStand/text
# @return location of subtext in text
# ----------------------------------------------------------------------
def _locateSubtext(subText, text):
    return [n for n in range(len(text)) if text.find(subText,
    n) == n]
# ----------------------------------------------------------------------

#Recursive Function
# ---peptide_combos-----------------------------------------------------
# Constructs all protein Combinations given a list of single rna peptide
# encode options
# @param reverse list of single rna peptide encodes
# @return list of combinations of protein
# ----------------------------------------------------------------------
def peptide_combos(peptide):
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
    print ('length of protein combinations with duplicates:',
    len(proteinCombinations))
    list(set(proteinCombinations))
    print ('length of unique protein combinations:', len(proteinCombinations),
    '\nlength of one protein:', len(proteinCombinations[0]))
    return proteinCombinations
# ----------------------------------------------------------------------

# ---_findAllSixRFs-----------------------------------------------------
# constructs list of all 6 reading frames of a dna text
# @param dna text
# #return list of all 6 reading frames
# ----------------------------------------------------------------------
def _findAllSixRFs(dna):
# first three reading frames with forward strand
    readingFrames = []
    readingFrames.append(dna)
    readingFrames.append(dna[1:])
    readingFrames.append(dna[2:])
# next 3 reading frames from reverse strand
    reverseComplement = reverse_complement(dna)
    readingFrames.append(reverseComplement)
    readingFrames.append(reverseComplement[1:])
    readingFrames.append(reverseComplement[2:])

    return readingFrames
# ----------------------------------------------------------------------

# ---translation--------------------------------------------------------
# Constructs an amino acid sequence by translating rna into amino acids
# @param rna string
# @return amino acid string
# ----------------------------------------------------------------------
def _translation(rna):
    codonList = [rna[i:i+3] for i in range(0, len(rna), 3)]

    protein = ''
    for codon in codonList:
        if codon in codonTable:
#        if codonTable[codon] == 'X':
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
    if type(dna) is list:
        for line in dna:
            reverseDna = _reverse(line)
# Applies complements of reversed input to dnaComplement
            for nucleotide in reverseDna:
                dnaComplement += _complement(nucleotide)
    else:
        reverseDna = _reverse(dna)
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
#byteText = urllib.request.urlopen(
#"https://stepik.org/media/attachments/lessons/96/RNA_codon_table_1.txt"
#).read()
# Convert bytes to string
#strText = (byteText.decode("utf-8"))
# Save text line by line as a list
#listText = strText.splitlines()

# Mark stop codons with an "X"
#for i in range(len(listText)):
#    if len(listText[i]) != 5:
#        listText[i] += 'X'

# Create Codon Table
codonTable = {"UUU" : "F", "UUC" : "F", "UUA" : "L", "UUG" : "L",
                  "UCU" : "S", "UCC" : "S", "UCA" : "S", "UCG" : "S",
                  "UAU" : "Y", "UAC" : "Y", "UAA" : "X", "UAG" : "X",
                  "UGU" : "C", "UGC" : "C", "UGA" : "X", "UGG" : "W",
                  "CUU" : "L", "CUC" : "L", "CUA" : "L", "CUG" : "L",
                  "CCU" : "P", "CCC" : "P", "CCA" : "P", "CCG" : "P",
                  "CAU" : "H", "CAC" : "H", "CAA" : "Q", "CAG" : "Q",
                  "CGU" : "R", "CGC" : "R", "CGA" : "R", "CGG" : "R",
                  "AUU" : "I", "AUC" : "I", "AUA" : "I", "AUG" : "M",
                  "ACU" : "T", "ACC" : "T", "ACA" : "T", "ACG" : "T",
                  "AAU" : "N", "AAC" : "N", "AAA" : "K", "AAG" : "K",
                  "AGU" : "S", "AGC" : "S", "AGA" : "R", "AGG" : "R",
                  "GUU" : "V", "GUC" : "V", "GUA" : "V", "GUG" : "V",
                  "GCU" : "A", "GCC" : "A", "GCA" : "A", "GCG" : "A",
                  "GAU" : "D", "GAC" : "D", "GAA" : "E", "GAG" : "E",
                  "GGU" : "G", "GGC" : "G", "GGA" : "G", "GGG" : "G"}
# Fill Codon Table with parsed information
#for i in range(len(listText)):
#    delimiter = ' '
#    codon = listText[i].split(delimiter)
#    for j in range(len(codon)):
#        delimiter = ' '
#        aminoAcid = codon[j].split(delimiter)
#        codonTable[codon[0]] = codon[1]

# Standard Input
rawData = sys.stdin.read().splitlines()
dna = rawData[0]
peptide = rawData[1]
myAnswer = (peptide_encode(dna, peptide))
#print (type(myAnswer))
print ('****My Answer:****')
for i in range(len(myAnswer)):
    item = myAnswer[i]
    if item[0] > 3:
        print ('\nreading frame:', item[0],
        '\npeptide:', item[1],
        '\nprotein:', _translation(_transcription(reverse_complement(item[1]))))
    else:
        print ('\nreading frame:', item[0],
        '\npeptide:', item[1],
        '\nprotein:', _translation(_transcription(item[1])))
answer = ['AAGGAAGTATTTGAGCCTCATTATTAC', 'AAAGAGGTGTTTGAACCTCATTACTAT',
'AAGGAGGTATTTGAACCCCACTATTAC',
'AAAGAAGTTTTCGAACCACATTATTAC',
'AAGGAAGTGTTTGAACCTCACTATTAT',
'AAAGAAGTTTTCGAGCCGCACTACTAC',
'AAGGAAGTATTCGAACCACATTACTAT',
'ATAATAATGCGGCTCGAATACTTCCTT',
'GTAGTAATGGGGCTCGAAAACCTCCTT',
'GTAGTAATGAGGTTCAAAAACCTCCTT',
'GTAGTAATGGGGTTCGAAGACTTCCTT',
'ATAATAGTGAGGCTCAAAAACTTCCTT',
'ATAGTAATGGGGTTCGAAGACTTCCTT',
'GTAGTAGTGCGGCTCAAAAACTTCCTT',
'ATAGTAATGAGGTTCGAAAACCTCTTT',
'ATAATAATGTGGCTCGAACACTTCTTT',
'GTAGTAATGGGGCTCAAACACCTCTTT',
'ATAGTAGTGAGGTTCGAAGACTTCCTT',
'GTAATAGTGCGGTTCAAAAACTTCCTT',
'ATAGTAGTGTGGTTCAAATACCTCCTT']
print ('\nAnswer:')
for peptide in answer:
    print ('peptide:', peptide, '\nprotein:', _translation(_transcription(peptide)))

answerList = []
for answer in myAnswer:
    answerList.append(answer[1])

for a in answerList:
    print (a)
