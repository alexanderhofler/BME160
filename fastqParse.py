#!/usr/bin/env python3
# Name: Your full name (CATS account username)
# Group Members: List full names (CATS usernames) or “None”

'''
Take user input in FASTQ form. Return the input on 7 different lines, each separated and categorized.
'''
import re
class FastqString (str):
    """
    Parses the single line of FASTQ text into several string outputs.
    """
    def _new_(self, sequence):
        ''' Returns the separated string without the @ symbol into the 7 different categories'''
        return str.__new__(self, sequence[1:])

    def parse(self):
        ''' Splits the string into 7 different strings, using the colon as the separating factor'''
        Seq = self.split (":")
        instr, rid, fcid, fcl, tnum, x, y = Seq
        print('Instrument = '+instr+'\nRun ID = '+rid+
              '\nFlow Cell ID = '+fcid+'\nFlow Cell Lane = '+fcl+
              '\nTile Number = '+tnum+'\nX-coord ='+x+'\nY-coord = '+y)

seq = input("Enter a single line of a FASTQ formatted file: ")
seq2 = FastqString(seq)
seq2.parse()
