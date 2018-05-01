#!/usr/bin/env python3
# Name: Alexander Hoefler (ahoefler)
# Group: none

class NucParams:
    """Creates 3 dictionaries to store fasta file information.
        Dict1: RNA codons
        Dict2: DNA codons
        Dict3: Nucleotides
    """

    rnaCodonTable = {
                    # RNA codon table
    # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }

    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    nucleotides = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N':0}

    def __init__(self, inString= ''):
        self.aaComp = {} #setting up empty dictionaries, values = 0
        self.codonComp = {}
        self.nucComp = {}

        for aa in ProteinParam.aa2mw:
            self.aaComp[aa] = 0

        for codon in self.rnaCodonTable:
            self.codonComp[codon] = 0

        for nuc in NucParams.nucleotides:
            self.nucComp[nuc] = 0

    def addSequence(self, inSeq):
        for nuc in inSeq:
            if nuc in NucParams.nucleotides.keys():  # Adds sequences using the inti method, checks that sequences entered are valid
                self.nucComp[nuc] += 1

        rnaSequence = inSeq.replace('T', 'U')  # Changes DNA sequence to RNA (T&U)

        for start in range(0, len(rnaSequence), 3):
            codon = rnaSequence[start: start + 3]
            if codon in self.rnaCodonTable:
                self.codonComp[codon] += 1  # Adds RNA sequence to dictionary w/ count.
                aa = self.rnaCodonTable[codon]
                if aa != '-':
                    self.aaComp[aa] += 1  # Adds amino acid to dictionary w/ count.


    def aaComposition(self):
        return self.aaComp

    def nucComposition(self):
        return self.nucComp #updates the acidComposition and returns it

    def codonComposition(self):
        return self.codonComp #Returns RNA with counts of codons

    def nucCount(self):
        nucTot = 0
        for count in self.nucComp.values():
            nucTot += count
        return nucTot #gives the total count of the nucleotides in the dictionary


class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):  # Recieved help from T.A. Logan during Lab
     """
     Takes the string splices it such that only valid characters are recieved.
     They are then read and the values are saved (not the string itself)
     """

    aaList = ''.join(protein).split()
    self.proteinInput = ''.join(aaList).upper()#return string as uppercase characters without invalid characters
    self.aaComp = {}

    for aa in self.aa2mw.keys():  # Finds and stores the values.
            self.aaComp[aa] = float(self.proteinInput.count(aa))

    def aaCount (self):
        """
        Counts and returns the number of amino acids in the sequence
        """

        aaTotal = 0
        for aa in self.proteinInput:
            if aa.upper() in self.aa2mw.keys():  # Checks if character in string is a valid aa.
                aaTotal += 1 #add 1 if the amino acid is part of the dictionary
            return aaTotal

    def pI (self):
        """
        Searches for which pH would make the sequence neutral (or close to)
        """

        loopCharge = 2**15  # Iteration invariant. Number is arbitrary but choose one that won't cut you short.
        bestPH = 0
        tempPH = 0
        while tempPH < 14.01: #checking that charge is within possible range (0-14)
            charge = abs(self._charge_(tempPH))
            if charge < loopCharge:
                loopCharge = charge
                best = tempPH
            tempPH += 0.01  # checking for best pH down to 2 decimals.
        return bestPH

    def aaComposition (self) :
        return self.aaComp #returns the dictionary setup at the beginning

    def _charge_ (self, pH):
        """
        Calculates the protein's specific pH's charge using n and c-terminus pKa's
        """

        posCharge = 0 #looking for total amount of positive charge
        for aa in self.aa2chargePos:
            numAminoAcid = self.aaComp[aa] #iterating through the dict to find the amount
            posCharge += numAminoAcid * ((10**self.aa2chargePos[aa])
            / (10**self.aa2chargePos[aa] + 10**pH))
        posCharge += (10**self.aaNterm) / (10**self.aaNterm + 10**pH)

        negCharge = 0 #same as with the positive charge but looking for negative charges
        for aa in self.aa2chargeNeg:
            numAminoAcids = self.aaComp[aa]
            negCharge += numAminoAcid * ((10**pH)
            / (10**self.aa2chargeNeg[aa] + 10**pH))
        negCharge += (10**pH) / (10**self.aaCterm + 10**pH)

        return posCharge - negCharge

    def molarExtinction (self):
        """
        We create a formula using amino acid absorbance at 280nm to find the molar extinction coefficient,
         which is used for the mass extinction coefficient
         """
        extinction = float(0)
        for aa in self.aa2abs280.keys():
            extinction += self.aaComp[aa] * self.aa2abs280[aa]

        return extinction


    def massExtinction (self):
        molecWeight =  self.molecularWeight()

        return self.molarExtinction() / molecWeight if molecWeight else 0.0

    def molecularWeight (self):
        """
        Calculates the proteins molecular weight after taking the loss of
        water from bond formation into account"""

        h2oWeight = 0
        water = self.mwH2O * (self.aaCount() - 1)
        for aa, count in self.aaComp.items():
            h2oWeight += (count * self.aa2mw[aa])  # Weight of the individual aa's - water released by bond formation
        return h2oWeight - water

# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))

        inString = input('protein sequence?: ')

if __name__ == "__main__":
    main()



import sys
class FastAreader :
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
