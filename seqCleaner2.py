#!/usr/bin/env python3
# Name: Alexander Hoefler (ahoefler)
# Group Members: none

class SequenceCleaner(str):
    """
    SequenceCleaner class returns a string object in upper case letters.
    Keyword arguments:
    sequence -- DNA sequence user enters
    """
    def __new__(self, sequence):
        """Returns a copy of sequence in upper case letters."""
        return str.__new__(self, sequence.upper())

    def getN(self):
        """Returns the number of N's in the sequence."""
        num_N = self.count("N")
        return '{'+str(num_N)+'}'

    def removeN(self):
        """Removes the ambiguous bases(denoted by “N”) from a sequence."""
        firstN = self.find("N")
        lastN  = self.rfind("N")
        print(self[:firstN]+self.getN()+self[lastN+1:])

"""
Builds a new string object of SequenceCleaner class based on what the user
inputs.
Removes the ambiguous parts of the sequence, outputs the “cleaned”
sequence, replacing the ambiguous parts with a count in {}’s
"""
DNA = input("Enter a DNA sequence: ")
seq = SequenceCleaner(DNA)
seq.removeN()
