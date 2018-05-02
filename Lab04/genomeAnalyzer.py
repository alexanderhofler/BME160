from sequenceAnalysis import NucParams,ProteinParam, FastAreader


class genomeAnalyzer:
    def __init__(self, filename='testGenome.fa'):
        """
        Sets the FastA file to be read
        """
        self.gc = self.gcContent(NucParams)
        self.genomeAnalysis(filename)
        self.myReader = FastAreader(filename)

    def genomeAnalysis(self, filename: object) -> object:
        """
        :param filename:
        :return: length and amino acid count
        """
        length, nucleotides = self.sequenceLength(self.myReader)
        print("sequence length = {:.2f}Mb".format(length))
        print(" ")
        print("GC content = {:.1f}%".format(self.gc))
        print(" ")

        codonCount = nucleotides.codonComposition() #count and sort codons
        codonCount.sort
        for codons in codonCount:
            total = codonCount[codons]
            aa = nucleotides.rnaCodonTable[codons]
            aaCount = nucleotides.aaComposition()[aa]

            if aaCount != 0:
                finalTotal = (total/aaCount) * 100
            else:
                finalTotal = (total/1) * 100

            print('{} : {} {:5.1f}% ({:6d})'.format(codons, aa, finalTotal, total))



    ##determining the GC Content
    def gcContent(self, nucleotides):

        nucComp = nucleotides.nucComposition()
        gc = nucComp['G'] + nucComp['C']
        gc = gc/nucleotides.nucCount() * 100

        return gc

    ##determining sequence length
    def sequenceLength(self, myReader):
        nucleotides = NucParams('')
        for head, seq in myReader.readFasta():
            nucleotides.addSequence(seq)

        length = nucleotides.nucCount()/1000000
        return length , nucleotides



genomeAnalyzer()
