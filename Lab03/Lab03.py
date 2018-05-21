
# coding: utf-8

# # Flow Control
# 
# 
#  - Due: April 23, 2018 at 11:55 pm, submit in Canvas
#  - program filename: proteinParams.py ( do not submit the Jupyter notebook)
#  - total points possible: 40
#  - extra-credit possible: 6

# ## Building and Testing
# For this lab, you will submit a single text file named proteinParams.py. You can do all of your work from within jupyter, then copy/paste your program to a new text file using your favorite text editor or python IDE. Test that program file !! On Mac systems you will use terminal, and on windows systems you will likely use cmd. You will then navigate to the directory where your proteinParam.py lives, then execute: <br>
# python3 proteinParam.py <br>
# or <br>
# python proteinParam.py <br>
# ## Protein parameters
# In this exercise, you’ll create a single Python program to calculate the physical-chemical properties of a protein sequence similar to what ProtParam outputs. Your task is to develop the ProteinParam class included in template form, below. The program includes the methods that you will need, along with the completed "main" that does all of the input and output.
# 
# When testing, you will type the protein sequence, hit return, the program will respond with the required output.  You can then enter a new protein string and a new analysis will be presented.  When you are finished, type ctrl-D to signal the end (this is done by holding the control key down as you type the letter D).
# 
# Your program will read in the protein sequence and print out the:
# 
#  - number of amino acids and total molecular weight,
#  - Molar extinction coefficient and Mass extinction coefficient,
#  - theoretical isoelectric point (pI), and
#  - amino acid composition
#  
# For example, if I enter the protein sequence:
# VLSPADKTNVKAAW 
# 
# then the program should output:

# <pre>
# Number of Amino Acids: 14
# Molecular Weight: 1499.7
# molar Extinction coefficient: 5500.00
# mass Extinction coefficient: 3.67
# Theoretical pI: 9.88
# Amino acid composition:
# A = 21.43%
# C = 0.00%
# D = 7.14%
# E = 0.00%
# F = 0.00%
# G = 0.00%
# H = 0.00%
# I = 0.00%
# K = 14.29%
# L = 7.14%
# M = 0.00%
# N = 7.14%
# P = 7.14%
# Q = 0.00%
# R = 0.00%
# S = 7.14%
# T = 7.14%
# V = 14.29%
# W = 7.14%
# Y = 0.00%
# </pre>

# # Hints:
# 
# The input sequence is not guaranteed to be uppercase and might contain unexpected characters.  
# Only count the following (A, C, D, E, F, G, H, I, L, K, M, N, P, Q, R, S, T, V, Y, W) or the lower-case equivalents, and ignore anything else. Math details are included below.
# 
# 
# To get full credit on this assignment, your code needs to:
# 
#  - Run properly (execute and produce the correct output)
#  - Include a docstring overview about what your program is designed to do with expected inputs and outputs
#  - include docstrings for every class, method within that class
#  - Include any assumptions or design decisions you made in writing your code as # comments
#  - Contain in-line comments using #-style where appropriate. Make sure to fix the template (below) to conform.
#  
#  - Submit your proteinParams.py file ( not the notebook) using Canvas.
# 
# Congratulations, you finished your third lab assignment!

# ## Design specification 
# 
# ### __init__
# 
# There are a number of ways to design this.  Your __init__ method could save an attribute which is just the input string. A more effective solution would compute and save the aaComposition dictionary. All of the protein parameter methods can operate very efficiently using a dictionary. Use of an aaComposition here will save you quite a bit of work.
# 
# ### aaCount()
# 
# This method will return a single integer count of valid amino acid characters found. Do not assume that this is the length of the input string, since you might have spaces or invalid characters that are required to be ignored.
# 
# ### aaComposition() - 4 points
# 
# This method is to return a dictionary keyed by single letter Amino acid code, and having associated values that are the counts of those amino acids in the sequence. Make sure to include all 20 amino acids.  Proper amino acids that are not represented in the sequence should have a value of zero. Note: if you have already calculated a composition dictionary in __init__, then just return that dictionary here.
# 
# ### molecularWeight() - 8 points
# 
# This method calculates the molecular weight (MW) of the protein sequence. If we have the composition of our protein, this is done by summing the weights of the individual Amino acids and excluding the waters that are released with peptide bond formation. 
# $$MW_{H_{2}O}+\sum^{aa}N_{aa}(MW_{aa}-MW_{H_{2}O})$$
# 
# ### \_charge\_(pH) -- 10 points
# 
# This method calculates the net charge on the protein at a specific pH (specified as a parameter of this method). The method is used by the pI method. I have marked it with the single _ notation to advise future users of the class that this is not part of the defined interface and it just might change.
# 
# If we have the composition of our protein, we can then calculate the net charge at a particular pH, using the pKa of each charged Amino acid and the Nterminus and Cterminus.
# 
# $$ netCharge = \left[\sum^{aa=(Arg,Lys,His,Nterm}N_{aa}\frac{10^{pKa(aa)}}{10^{pKa(aa)}+10^{pH}}\right] - 
#                \left[\sum^{aa=(Asp,Glu,Cys,Tyr,Cterm}N_{aa}\frac{10^{pH}}{10^{pKa(aa)}+10^{pH}}\right]
# $$
# 
# I have provided pKa tables for each AA, and the pKa for the N-terminus and C-terminus.
# 
# ### pI() - 10 points
# 
# The theoretical isolelectric point can be estimated by finding the  particular pH that yields a neutral net Charge (close to 0). There are a few ways of doing this, but the simplest might be to iterate over all pH values to find the one that is closest to 0. Doing this by hand is painful, but its not that bad to do computationally.  Remember that we want to find the  best pH, accurate to 2 decimal places.
# 
# #### extra credit (3 points) for pI() method
# 
# Another way of doing the pI calculation would use a binary search over the pH range. This works because we expect a single zero crossing to exist in the range, and the function will be well behaved across the range of charge() as a function of pH (0-14 range). You then make the algorithm operate to any specified precision using an optional parameter (see Model page 33 for an example - set the default parameter: precision = 2).
# 
# ### molarExtinction() - 8 points
# 
# The extinction coefficient indicates how much light a protein absorbs at a certain wavelength. It is useful to have an estimation of this coefficient for measuring a protein with a spectrophotometer at a wavelength of 280nm. It has been shown by Gill and von Hippel that it is possible to estimate the molar extinction coefficient of a protein from knowledge of its amino acid composition alone. From the molar extinction coefficient of tyrosine, tryptophan and cystine at a given wavelength, the extinction coefficient of the native protein in water can be computed using the following equation.
# 
# $$E = N_Y E_Y + N_W E_W + N_C E_C$$
# where:
# 
# - $N_Y$ is the number of tyrosines, $N_W$ is the number of tryptophans, $N_C$ is the number of cysteines,
# - $E_Y$ , $E_W$, $E_C$ are the extinction coefficients for tyrosine, tryptophan, and cysteine respectively.
# 
# I have supplied the molar extinction coefficients at 280nm for each of these residues in a dictionary(aa2abs280) in the program template.
# 
# Note that we will assume for this exercise that all Cysteine residues are represented as Cystine. Under reducing conditions, Cystine does not form however and Cysteine residues do not contribute to absorbance at 280nm.
# 
# ### massExtinction()
# 
# We can calculate the Mass extinction coefficient from the Molar Extinction coefficient by dividing by the molecularWeight of the corresponding protein.
# 
# #### extra credit for molarExtinction() and massExtinction() - 3 points
# 
# As mentioned above, we are assuming that all Cysteine residues are present as Cystine.  Provide an optional parameter to both molarExtinction() and massExtinction() to calculate these values assuming reducing conditions. Use an optional parameter with a default of True (Cystine=True) to allow evaluation of both molar and mass extinction under both oxidizing (default) and reducing conditions.  (see Model page 33)

# ## Protein Param

# In[ ]:


#!/usr/bin/env python3
# Name: Alexander Hoefler (ahoefler)
# Group Members: None

"""
Program meant to calculate the properties of the amino acid sequence inputted by the user
"""
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
"""
sets up dicitonary and defines valid aa's as well as returning
the inputted string with only uppercase characters"""
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

    def __init__ (self, protein):  # Recieved help from T.A. Logan during Lab
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

    def molarExtinction (self, Cysteine = True):
       """
       Provides the exctiction coefficient using tyrosine, tryptophans, and cysteines.
       Does so at a wavelength of 280nm (each provided in the dictionary above)
       """
        tyro = self.aaComp['Y'] * self.aa2abs280['Y']
        tryp = self.aaComp['W'] * self.aa2abs280['W']
        cyst = self.aaComp['C'] * self.aa2abs280['C']
        extinction = tyro + tryp + cyst
        return extinction

    def massExtinction (self, Cysteine = True):
        molecWeight = self.molecularWeight()
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


# In[7]:


VLSPADKTNVKAAW


# In[ ]:


VLSPADKTNVKAAW

