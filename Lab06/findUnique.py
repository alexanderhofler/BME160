#!/usr/bin/env python3
# Alexander Hoefler (ahoefler)
# Group: none

"""
    Input for STDIN (in order to compile):
    python3 findUnique.py <[bos-tRNA.fa].fa (or whatever the name of your preferred file is)
"""
import sequenceAnalysis as seqAnal


class findUnique:
    def __init__(self):
        """
        Reads tRNA from the inputted Fasta file. Adds the power set of each sequence to a list

        """
        self.powerSetList = []
        self.uniqueList = []  # List fo tRNA substrings
        self.headerSequenceDictionary = {}  # sequence saving dict.
        fastaFile = seqAnal.FastAreader()  # pulls in fasta reader and initializes here to read in fasta file
        count = 0

        for header, sequence in fastaFile.readFasta():
            filteredSequence = self.removeChar(sequence)  # makes sure only valid characters are being used
            self.headerSequenceDictionary[count] = [header, filteredSequence]
            mypowerSet = self.powerSet(filteredSequence)
            self.powerSetList.append(mypowerSet)  # Add power set (after finding) for each tRNA to a list.
            count += 1

    def removeChar(self, sequence):
        """
        Removes unwanted characters from the sequence
        """
        noDashSequence = sequence.replace('-','')
        noUnderscoreSequence = noDashSequence.replace('_','')
        return noUnderscoreSequence

    def powerSet(self, sequence):
        """
        Returns tRNA substrings
        """
        powerSets = set()
        for index in range(len(sequence)):
            size = len(sequence)
            while size > index:
                powerSets.add(sequence[index:size])
                size -= 1
        return powerSets #goes through the sequence piecewise before returning power values

    def findUniques(self):
        """
        Looks for duplicate tRNA sequences and removes from all but the appriopiate set to
        obtain a unique set of tRNA subsequences.
        """
        for allSet in self.powerSetList:
            union = set()
            copySet = set()
            copyList = self.powerSetList.copy()
            copySet = allSet.copy()  # Make duplicate of the possible pwer set
            copyList.remove(copySet)  # Remove power set from list.
            for powerSet in copyList: # removes the duplicate elements between union and copySet
                union = union.union(powerSet)  # The union of all the other tRNA sets.
            copySet.difference_update(union)  # Update the copySet, removing elements found in union.

            newSet = copySet.copy()
            for string1 in copySet:  # checks prior strings to identify more possible substrings
                uniqueSet = copySet.copy()
                uniqueSet.remove(string1)
                for string2 in uniqueSet:
                    if string1 in string2:
                        if len(string1) < len(string2):
                            newSet.discard(string2)  # Gets rid of teh larger substrings.

            self.uniqueList.append(newSet)

    def printSequences(self):
        """
        Prints to STDOUT commandline
        """
        for index in range(0,len(self.headerSequenceDictionary)):
            mySequence = self.headerSequenceDictionary[index]
            header = mySequence[0]
            sequence = mySequence[1]
            print(header)
            print(sequence)
            size = len(sequence)
            for position in range(0, size):  # Note positions of the substrings for later
                for substring in self.uniqueList[index]:
                    substringLength = len(substring)
                    if substring == sequence[position:position + substringLength]:
                        output = ('.')*position + substring  # ipnut appropiate amount of "."'s then add found substring
                        print(output)

def main():
    """
    Runs above code to find and print header, sequence and substrings from the fasta file.
    """
    mytRNA = findUnique()
    mytRNA.findUniques()
    mytRNA.printSequences()

if __name__ == "__main__":
    main()
