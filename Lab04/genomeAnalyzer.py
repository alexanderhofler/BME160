#!/usr/bin/env python
"""Analyze Genome """
########################################################################
# File: genomeAnalyzer.py
#  executable: genomeAnalyzer.py

# Author: Alex Hofler
# History: 06/05/17 Created
########################################################################

from sequenceAnalysis import NucParams, FastAreader, ProteinParam


class genomeAnalyzer:
    # TODO you need docstrings here like I did for sequenceAnalysis
    def __init__(self, genome_path='testGenome.fa'):
        #TODO need docstrings
        self.genome = genome_path
        self.myReader = FastAreader(self.genome)
        self.nuc_params = NucParams()
        for head, seq in self.myReader.readFasta():
            self.nuc_params.addSequence(seq)

    def genomeAnalysis(self):
        """Compute and print sequence length, gc content and amino acid composition"""
        length = self.sequenceLength()
        # Alex, you should not change the length of the sequence until you format it because if you need to use the
        # function for something else you would have to change a bunch of stuff for it to work

        print("sequence length = {:.2f}Mb\n".format(length/1000000.0))
        # same idea goes for gcContent. It should return the gcContent and if you want to convert to percentage, you can
        gc = self.gcContent()
        print('GC content = {:.1f}%\n'.format(gc*100))

        codon_counts = self.nuc_params.codonComposition()
        aa_comp = self.nuc_params.aaComposition()

        # go through each amino acid
        for aa in sorted(aa_comp):
            # created another data structure aaTable in nuc_params to deal with going from aa to RNA codon
            for codon in sorted(self.nuc_params.aaRnaTable[aa]):
                codon_total = codon_counts[codon]
                aa_count = aa_comp[aa]
                # calculate relative codon usage for each codon and print
                if aa_count != 0:
                    total = (codon_total/aa_count) * 100
                else:
                    print("This happens")
                    total = (codon_total/1) * 100

                print('{} : {} {:5.1f}% ({:6d})'.format(codon, aa, total, codon_total))

    def gcContent(self):
        """Compute GC content of entire fasta file"""
        nuc_comp = self.nuc_params.nucComposition()
        gc = nuc_comp['G'] + nuc_comp['C']
        gc = gc/self.nuc_params.nucCount()
        return gc

    def sequenceLength(self):
        """Computes the sequence length of a given genome
        :return length, nucParams class object
        """
        return self.nuc_params.nucCount()


# TODO you should have a main function here and pass in the argument
def main():
    """Analyzes genome"""
    genome = 'testGenome.fa'
    gA = genomeAnalyzer(genome_path=genome)
    gA.genomeAnalysis()

if __name__ == '__main__':
    main()
