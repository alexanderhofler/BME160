#!/usr/bin/env python3
# Name: Alexander Hoefler (ahoefler)
# Group: None - Assistance by Noah Dove (tutor)
#
# ORFfinder class created in sequenceAnalysis:
#    includes start and stop positioning as well as the extra credit options
#       --> slight edit made to the argparse provided in order to
#                  set default minGene (mG to 100)
# init
#     read in sequence, and create list of ORFs found
#
# findORF: look for codons stepwise, slice to look every 3 starting from frame 0
#           if start codon is found
#               save the position under appropiate list
#           if stop codon found before/without start codon
#               assume start = 1
#           if start and stop codons are found
#               start = beginning of the list [0]
#               stop = beginning of stop codon (i) + 3
#           if a start codon is found without a stop codon
#               start = beginning of the list [0],
#               assume entire length of the sequence involved
#
# for the reverse strand: do the same but in the reverse direction (look at the compliment strand)
#
# saveOrf: add the start, stop, length, and the frame of the ORF to the list
#
# reverseComplement: create a reverse compliment of the sequence
#     helper method to take the reverse complement of DNA seq.
#
# main method:
# read in sequence from fasta file,
#   prints frame, start, stop, length, and takes input from commandline on variations

from sequenceAnalysis import OrfFinder, FastAreader

class CommandLine():

    """
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was

    set as an option using add_argument, then myCommandLine.args.requiredbool will

    name that option.
    """
    def __init__(self, inOpts=None):
        """
        Implement a parser to interpret the command line argv string using argparse.
            --> Argparse adjusted so that default is minGene = 100
        """
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Finds the largest ORF in a DNA sequence',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store_true', help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(100, 200, 300, 500, 1000), action='store',
                                 help='minimum Gene length', default = 100)
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


########################################################################
# Main
# Here is the main program
#
#
########################################################################

def main(inCL=None):
    """Reads in a fasta file and outputs the ORFs frame, start, stop, and length position on a output file."""
    if inCL is None:
        myCommandLine = CommandLine()
        if myCommandLine.args.longestGene:
            fastaFile = FastAreader()
            for header, sequence in fastaFile.readFasta():
                print(header)
                orfData = OrfFinder(sequence)
                orfData.findOrfs()
                orfData.findRevOrfs()
                filteredList = filter(lambda orf: orf[3] > myCommandLine.args.minGene, orfData.orfs)  # Filters out the ORFS depending on the minGene arg.
                for frame, start, stop, length in sorted(filteredList, key=lambda orf: orf[3], reverse = True):  # Sorts the list of ORFs by length.
                    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length))
    else:
        myCommandLine = CommandLine(inCL)
    print(myCommandLine.args)

if __name__ == "__main__":
    main()
