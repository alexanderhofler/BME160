#!/usr/bin/env python3
# Name: Alexander Hoefler (ahoefler)
# Group: None - Assistance from Noah Dove (tutor)

# Adjusted/corrected for lab05, ORF class added
import sys

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

    validNucleotides = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N':0}

    def __init__(self):
        """
        Sets up dictionaries with keys, and sets values to zero.
        """
        self.aminoAcidComposition = {}
        self.codonsComposition = {}
        self.nucleotideComposition = {}

        for aa in ProteinParam.aa2mw:  # Initializes Amino Acid dictionary.
            self.aminoAcidComposition[aa] = 0

        for codon in self.rnaCodonTable:  # Initializes Codon Composition dictionary.
            self.codonsComposition[codon] = 0

        for nuc in NucParams.validNucleotides:  # Initializes Nucleotide Composition dictionary.
            self.nucleotideComposition[nuc] = 0

    def addSequence(self, thisSequence):
        """
        Adds all data (new and old) to the init method
        """
        for nuc in thisSequence:
            if nuc in NucParams.validNucleotides.keys():  # Adds sequences using the inti method, checks that sequences entered are valid
                self.nucleotideComposition[nuc] += 1
        rnaSequence = thisSequence.replace('T', 'U')  # Changes DNA sequence to RNA (T&U)

        for start in range(0, len(rnaSequence), 3):
            codon = rnaSequence[start: start + 3]
            if codon in self.rnaCodonTable:
                self.codonsComposition[codon] += 1  # Adds RNA sequence to dictionary w/ count.
                aa = self.rnaCodonTable[codon]
                if aa != '-':
                    self.aminoAcidComposition[aa] += 1  # Adds amino acid to dictionary w/ count.


    def aaComposition(self):
        return self.aminoAcidComposition

    def nucComposition(self):
        return self.nucleotideComposition #updates the acidComposition and returns it

    def codonComposition(self):
        return self.codonsComposition #Returns RNA with counts of codons

    def nucCount(self):
        nucTotal = 0
        for count in self.nucleotideComposition.values():
            nucTotal += count
        return nucTotal #gives the total count of the nucleotides in the dictionary


class ProteinParam(str):
    """Creates a ProteinParam string object in upper case letters.
    Has several functions that check work to perform chemical and
    physical property analysis.
    """

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):
        """Takes in sequence from users input then creates a list of
        strings and concatnates to create a string object.
        Then creates aaDictionary with the keys being the valid aa's of
        users input and the value of each one being their corresponding count.
        Args:
            param1 (str): String of amino acids.
        """
        myList = ''.join(protein).split()
        self.proteinString = ''.join(myList).upper()
        self.aaDictionary = {}

        for aa in self.aa2mw.keys():  # Iterates through aa2mw dictionary, the valid aa's.
            self.aaDictionary[aa] = float(self.proteinString.count(aa))  # Stores values.

    def aaCount(self):
        """
        Counts and returns the number of amino acids in the sequence.
        Ignores invalid aa's,
        """
        aaTotal = 0
        for aa in self.proteinString:
            if aa.upper() in self.aa2mw.keys():  # Checks if character in string is a valid amino acid.
                aaTotal += 1
        return aaTotal

    def pI(self):
        """
        Searches for which pH would make the sequence neutral (or close to).
        --> Looks up to 2 decimal points.
        Returns:
            (float): Best pH.
        """
        bigCharge = 2 ** 15
        bestPH = 0
        particularPH = 0
        while particularPH < 14.01:
            charge = abs(self.__charge__(particularPH))  # 0-14 (normal pH range)
            if charge < bigCharge:
                bigCharge = charge
                bestPH = particularPH
            particularPH += 0.01  # Iterates to 2 decimal points for pH.
        return bestPH

    def aaComposition(self):
        """
        Returns:
            aaDictionary created in __init__ method.
        """
        return self.aaDictionary #returns the dictionary setup at the beginning

    def __charge__(self, pH):
        """
        Calculates the protein's specific pH's charge using n and c-terminus pKa's
        """

        posCharge = 0
        for aa in self.aa2chargePos:  #looking for total amount of positive charge
            nAA = self.aaDictionary[aa]
            posCharge += nAA * ((10 ** self.aa2chargePos[aa]) #iterating through the dict to find the amount
                                / (10 ** self.aa2chargePos[aa] + 10 ** pH))
        posCharge += (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)

        negCharge = 0 #same as with the positive charge but looking for negative charges
        for aa in self.aa2chargeNeg:
            nAA = self.aaDictionary[aa]
            negCharge += nAA * ((10 ** pH)
                                / (10 ** self.aa2chargeNeg[aa] + 10 ** pH))
        negCharge += (10 ** pH) / (10 ** self.aaCterm + 10 ** pH)

        netCharge = posCharge - negCharge

        return netCharge

    def molarExtinction(self):
        """Create a formula using amino acid absorbance at 280nm to find the molar extinction coefficient,
         which is used for the mass extinction coefficient
        Returns:
            (float): Molar extinction coefficient.
        """
        tyrosine = self.aaDictionary['Y'] * self.aa2abs280['Y']
        tryptophans = self.aaDictionary['W'] * self.aa2abs280['W']
        cysteines = self.aaDictionary['C'] * self.aa2abs280['C']
        molarEx = tyrosine + tryptophans + cysteines
        return molarEx

    def massExtinction(self):
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        """
        Calculates the proteins molecular weight after taking the loss of
        water from bond formation into account
        Returns:
            (float): Molecular weight.
        """
        aaWeight = 0
        h2o = self.mwH2O * (self.aaCount() - 1)
        for aa, count in self.aaDictionary.items():
            aaWeight += (count * self.aa2mw[aa])   # Weight of the individual aa's - water released by bond formation
        return aaWeight - h2o


class FastAreader:

    def __init__ (self, fileName=''):
        '''contructor: saves attribute fileName '''
        self.fileName = fileName

    def doOpen (self):
        if self.fileName is '':
            return sys.stdin
        else:
            return open(self.fileName)

    def readFasta (self):

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

class OrfFinder():
    """
    Takes in a sequence of a FASTA file and stores ORFs in a lists of lists.
    Attributes:
        attr1 (list): List of valid stop codons.
        attr2 (list): List of valid start codon.
        attr3 (dict): Dictionary of valid nucleotides.
    """
    stop_codons = ['TGA', 'TAG', 'TAA']
    start_codons = ['ATG']
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, seq):
        """Sets up list containing many lists of oRFs found in the fasta sequence.
        Contains the frame, start, stop, length of the ORF.
        Args:
            param1 (str): String of DNA from FASTA file.
        """
        self.seq = seq
        self.orfs = []  # The lists where the ORFs frame, start, stop, and length are sent/setup/

    def findOrfs(self):
        """
        Find Orfs on 3'-5' strand and return list of Orfs.
        """
        start_positions = []
        foundStart = 0
        foundCodon = 0

        for frame in range(0,3):  # Check first 3 frames
            foundStart = 0  # After finding codons and start codons, notify.
            foundCodon = 0
            start_positions = []  # Clears the start position list for each frame.
            for i in range(frame, len(self.seq), 3):
                codon = self.seq[i:i+3] # Sets length of the codon to 3 nucs.
                if codon == 'ATG':  # default start codon - can be adjusted in the command line.
                    start_positions.append(i)
                    foundStart = 1
                    foundCodon = 1

                if codon in OrfFinder.stop_codons and foundStart:
                    start = start_positions[0] + 1 - frame
                    stop = i + 3
                    length = stop - start + 1
                    self.saveOrf((frame%3) + 1, start, stop, length)
                    start_positions = []
                    foundStart = 0
                    foundCodon = 1

                if not foundCodon and codon in OrfFinder.stop_codons:  # If no start codon was found but stop codon found. - dangling stop
                    start = 1
                    stop = i + 3
                    length = stop - start + 1
                    self.saveOrf((frame%3) + 1, start, stop, length)
                    start_positions = []
                    foundCodon = 1

            if foundStart:  # If no stop codon was found but start codon was found. - dangling start
                start = start_positions[0] + 1
                stop = len(self.seq)
                length = stop - start + 1
                self.saveOrf((frame%3) + 1, start, stop, length)

        return self.orfs

    def findRevOrfs(self): #go trhough test file and comment out working lines, then compare orfs by hand to find error
        """ Find Orfs on 5'-3' strand and returns that list of Orfs.
        """
        reverseComp = self.reverseComplement()
        start_positions = []
        foundStart = 0
        foundCodon = 0

        for frame in range(0, 3):  # Check  frames 1, 2, 3
            foundStart = 0  # Flag when finding codons and start codons.
            foundCodon = 0
            start_positions = []  # Clears the start position list for each frame.
            for i in range(frame, len(reverseComp), 3):
                codon = reverseComp[i:i + 3]  # The codon is 3 nucleotides.

                if codon == 'ATG':  # When start codon is found.
                    start_positions.append(i)
                    foundStart = 1
                    foundCodon = 1

                if codon in OrfFinder.stop_codons and foundStart:  # We have a start and stop codon.
                    stop = len(reverseComp) - start_positions[0]
                    start = len(reverseComp) - (i+2)
                    if frame == 1: stop += 1
                    elif frame == 2: stop += 2
                    length = stop - start + 1
                    self.saveOrf(-1 * ((frame%3) + 1), start, stop, length)
                    start_positions = []
                    foundStart = 0
                    foundCodon = 1

                if not foundCodon and codon in OrfFinder.stop_codons:  # If no start codon was found but stop codon found. - dangling stop (again)
                    start = len(reverseComp) - i - 2
                    stop = len(reverseComp)
                    length = stop - start + 1
                    self.saveOrf(-1 * ((frame%3) + 1), start, stop, length)
                    start_positions = []
                    foundCodon = 1

            if foundStart:  # If no stop codon was found but start codon was found. - dangling start (again)
                start =  start_positions[0] + 1
                stop = 1
                length = stop - start + 1
                self.saveOrf(-1 * ((frame%3) + 1), start, stop, length)

        return self.orfs

    def saveOrf(self, frame, start, stop, length):
        """ Saves ORF info
        """
        self.orfs.append([frame, start, stop, length])

    def reverseComplement(self):
        """ Returns the reversed and complimentary strand of DNA
        """
        return ''.join([self.complement[base] for base in self.seq[::-1]])  # Dictionary to find reverse complement of DNA sequence.
