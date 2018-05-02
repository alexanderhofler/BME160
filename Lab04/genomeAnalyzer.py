from sequenceAnalysis import NucParams, FastAreader, ProteinParam


class genomeAnalyzer:
    def __init__(self, filename='testGenome.fa'):
        self.genomeAnalysis(filename)

    def genomeAnalysis(self, filename: object) -> object:

        self.myReader = FastAreader(filename)

        length, nucParams = self.sequenceLength(self.myReader)
        print("sequence length = {:.2f}Mb".format(length))
        print(" ")

        self.gc = self.gcContent(nucParams)

        print('GC content = {:.1f}%'.format(self.gc))
        print(" ")

        codonCount = nucParams.codonComposition()

        for codons in sorted(codonCount): # sort codons in alpha order, by amino Acid
            codonTotal = codonCount[codons]
            aa = nucParams.rnaCodonTable[codons]
            aaCount = nucParams.aaComposition()[aa]
            # calculate relative codon usage for each codon and print
            if aaCount != 0:
                total = (codonTotal/aaCount) * 100
            else:
                total = (codonTotal/1) * 100

            print('{} : {} {:5.1f}% ({:6d})'.format(codons, aa, total, codonTotal))

    def gcContent(self, nucParams):
        nucComp = nucParams.nucComposition()
        gc = nucComp['G'] + nucComp['C']
        gc = gc/nucParams.nucCount() * 100

        return gc

    def sequenceLength(self, myReader):
        nucParams = NucParams('')
        for head, seq in myReader.readFasta():
            nucParams.addSequence(seq)

        length = nucParams.nucCount()/1000
        return length , nucParams



genomeAnalyzer()
