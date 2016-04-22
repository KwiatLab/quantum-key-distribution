# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 2012
"""
#from pylab import *
import numpy
def log2(x):
    return numpy.log2(x)

"""
Given a sequence with an alphabet size 'alphabet' with each letter equally likely (symmetric), the entropy per letter is:
1) entropy_per_letter =  -log_2(1/alphabet)

Given probability of having a coincidence in 2 sequences as 'p1g1', then the probability of the second sequence
having the same letter is p1g1.
The probability of the second sequence having a specific different letter is pdifferent = (1-p1g1)/(alphabet-1).

This gives the entropy of sequence 2 given sequence 1 as:
2) entropy_s2g1 =   -p1g1*log2(p1g1)-(alphabet-1)*pdifferent*log2(pdifferent)

Subtracting the two gives the maximum possible entropy left over after running a nonbinary SW code on this sequence
"""

def theoretical_nb(p1g1,alphabet):
    #p1g1 is coincidence
    #p0g1 is failcoincidence
    print "pCoincidence:",p1g1
    print "Alphabet:",alphabet
    
    #Each letter of the alphabet has the same probability
    p0g1=(1-p1g1)/(alphabet-1)
    p_letter = 1/float(alphabet)
    
    entropy_per_letter = -log2(p_letter)
    print "\nEntropy Per Letter:",entropy_per_letter
    
    entropy_s2gs1 = -p1g1*log2(p1g1)-(alphabet-1)*p0g1*log2(p0g1)
    entropy_left = entropy_per_letter-entropy_s2gs1
    
    print "Minimum Entropy Sent:",entropy_s2gs1
    print "Maximum Entropy Left:",entropy_left
    
    return entropy_left

"""
Given a sequence with an alphabet size 'alphabet' with each letter equally likely (symmetric), the entropy per letter is:
1) entropy_per_letter =  -log_2(1/alphabet)

Given probability of having a coincidence in 2 sequences as 'p1g1', then the probability of the second sequence
having the same letter is p1g1.
The probability of the second sequence having a specific different letter is pdifferent = (1-p1g1)/(alphabet-1).

This gives the entropy of sequence 2 given sequence 1 as:
2) entropy_s2g1 =   -p1g1*log2(p1g1)-(alphabet-1)*pdifferent*log2(pdifferent)

Subtracting the two gives the maximum possible entropy left over after running a nonbinary SW code on this sequence
"""


def theoretical_mat_nb(p_letter,transmat,alphabet):
    #p1g1 is coincidence
    #p0g1 is failcoincidence
    print "Transition Matrix:",transmat
    print "Alphabet:",alphabet
    
    #Each letter of the alphabet has the same probability
    #p_letter = 1/float(alphabet)
    
    p_letA = numpy.dot(p_letter,transmat.transpose())
    print("PA:",p_letA)
    print("PB:",p_letter)
    entropy_per_letter = 0.0
    for i in xrange(len(p_letA)):
        if (p_letA[i] > 0):
            entropy_per_letter-=p_letA[i]*log2(p_letA[i])
    print "\nEntropy Per Letter:",entropy_per_letter
    
    entropy_s2gs1 = 0.0
    for i in xrange(len(transmat[:,0])):
        for j in xrange(len(transmat[0,:])):
            if (transmat[i,j] > 0.0 and p_letter[j] > 0.0):
                entropy_s2gs1 -= p_letter[j]*transmat[i,j]*log2(transmat[i,j])
    #entropy_s2gs1/=float(len(transmat[0,:]))
    
    entropy_left = entropy_per_letter-entropy_s2gs1
    
    if (entropy_left<=0.0):
        entropy_left =0.0

    print "Minimum Entropy Sent:",entropy_s2gs1
    print "Maximum Entropy Left:",entropy_left
    
    return entropy_left

"""
This calculates the same thing as the above function - except for a binary code

Given p1 as the probability of a 1 in a binary sequence
p0 = 1-p1

1) entropy_per_bit = -p1*log2(p1)-p0*log2(p0)

Given p1g1 as the probability of 1 in another sequence given this sequence (assuming that both sequences have the same statistics):
p0g1 = 1-p1g1
p1g0=p0g1*p1/p0
p0g0=1-p1g0

Then the entropy of the unknown sequence given a 1 in this sequence is
entropy_g1=-p1g1*log2(p1g1)-p0g1*log2(p0g1)

And entropy of unknown sequence given a 0 in this sequence is
entropy_g0=-p1g0*log2(p1g0)-p0g0*log2(p0g0)

Giving a total entropy per bit in the unknown sequence given this sequence:
2) entropy_s2gs1 = entropy_g1*p1+entropy_g0*p0

Subtracting the two gives the maximum possible entropy left over after running a binary SW code on this sequence

"""
def theoretical(p1a,p1g1a,p1b,p1g1b):
    
    #Probability calculations
    p0a=1-p1a
    p0b=1-p1b
    p0g1a = 1-p1g1a
    p0g1b = 1-p1g1b
    
    p1g0a=p0g1a*p1a/p0a
    p1g0b=p0g1b*p1b/p0b
    p0g0a=1-p1g0a
    p0g0b=1-p1g0b
    
    print "A: p0   %f\tp1   %f\np0g1 %f\tp1g1 %f\np0g0 %f\tp1g0 %f\n"%(p0a,p1a,p0g1a,p1g1a,p0g0a,p1g0a)
    print "B: p0   %f\tp1   %f\np0g1 %f\tp1g1 %f\np0g0 %f\tp1g0 %f\n"%(p0b,p1b,p0g1b,p1g1b,p0g0b,p1g0b)
    
    #Entropy calculations
    entropy_per_bit_a = -p1a*log2(p1a)-p0a*log2(p0a)
    entropy_per_bit_b = -p1b*log2(p1b)-p0b*log2(p0b)
    
    print "A: Entropy Per Bit:",entropy_per_bit_a
    print "B: Entropy Per Bit:",entropy_per_bit_b
    
    entropy_g1a=-p1g1a*log2(p1g1a)-p0g1a*log2(p0g1a)
    entropy_g0a=-p1g0a*log2(p1g0a)-p0g0a*log2(p0g0a)
    
    entropy_agb = entropy_g1a*p1a+entropy_g0a*p0a
    
    entropy_left = entropy_per_bit_a-entropy_agb
    if (entropy_left<=0.0):
        entropy_left =0.0

    print "Minimum Entropy Sent:",entropy_agb
    print "Maximum Entropy Left:",entropy_left
    
    if (numpy.isnan(entropy_left)): return 0.0
    return entropy_left



def entropy_calculate(
    p1a=0.05,        #The probability of a photon in a time bin (p1)
    p1g1a=0.4,       #The data's heralding efficiency (coincidence rate) (p1 given 1)
    p1b=0.05,        #The probability of a photon in a time bin (p1)
    p1g1b=0.4,       #The data's heralding efficiency (coincidence rate) (p1 given 1)
    alph=16,        #The alphabet (frame size used)
    b_p1a=0.05,      #Probability of 1 in binary code
    b_p1g1a=0.85,    #Coincidence rate of binary slepian-wolf 1s
    b_p1b=0.05,      #Probability of 1 in binary code
    b_p1g1b=0.85,    #Coincidence rate of binary slepian-wolf 1s
    nb_p1g1=0.7,    #Probability of coincidence in nonbinary code
    nb_bperf = 64,  #Effective number of original sequence bits per nonbinary letter
    b_test = 1.0,   #Percentage sent in parities
    nb_test = 0.8,  #Percentage sent in parities nonbinary
    nb_mat=None     #Transition matrix for nonbinary code
    ):

    print "\n\nTHEORY:\n\n"

    print "TOTAL SEQUENCE THEORETICAL ENTROPY:"
    ideal_entropy = theoretical(p1a,p1g1a,p1b,p1g1b)

    print "\n\nSLEPIAN-WOLF THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is per alphabet bits of the original sequence
    b_entropy = theoretical(b_p1a,b_p1g1a,b_p1b,b_p1g1b)/alph

    print "\nPercentage of total entropy recovered:",b_entropy/ideal_entropy

    print "\n\nNON-BINARY THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is for a large amount of original sequence bits
    nb_entropy = 0.0
    
    if (nb_mat!=None):
        nb_entropy = theoretical_mat_nb(nb_mat,alph)/nb_bperf
    else:
        nb_entropy = theoretical_nb(nb_p1g1,alph)/nb_bperf

    print "\nPercentage of total entropy recovered:",nb_entropy/ideal_entropy

    #These are calculations of entropy retained for the actual codes
    test_b_entropy = (1-b_test)/alph
    test_nb_entropy = log2(alph)*(1-nb_test)/nb_bperf

    print "\n\nCURRENT TEST RESULTS:"
    print "              Entropy   P_theory P_total"
    print "Slepian-Wolf: %f  %f %f"%(test_b_entropy,test_b_entropy/b_entropy,test_b_entropy/ideal_entropy)
    print "Non-Binary:   %f  %f %f"%(test_nb_entropy,test_nb_entropy/nb_entropy,test_nb_entropy/ideal_entropy)
        
    print "\n\nTOTAL PERCENTAGE OF ENTROPY RETAINED:"
    print "Theory: %f"%((nb_entropy+b_entropy)/ideal_entropy)
    print "Test:   %f"%((test_nb_entropy+test_b_entropy)/ideal_entropy)

    print "Test/Theory",(test_nb_entropy+test_b_entropy)/(nb_entropy+b_entropy)
    
    print "\n\n"
    
    return (ideal_entropy,b_entropy,nb_entropy)

def entropy_calculate2(
    p1a=0.05,        #The probability of a photon in a time bin (p1)
    p1g1a=0.4,       #The data's heralding efficiency (coincidence rate) (p1 given 1)
    p1b=0.05,        #The probability of a photon in a time bin (p1)
    p1g1b=0.4,       #The data's heralding efficiency (coincidence rate) (p1 given 1)
    alph=16,        #The alphabet (frame size used)
    b_plet=None,    #Probability of each letter
    b_mat=None,     #Transition matrix for binary code
    nb_bperf = 64,  #Effective number of original sequence bits per nonbinary letter
    nb_plet=None,   #nonbinary probability of each letter
    nb_mat=None     #Transition matrix for nonbinary code
    ):

    print "\n\nTHEORY:\n\n"

    """
    The following are based upon actual sequences
    """


    #p1_bsw = 1-(1-p1_data)**alphabet    #The probability of having at least one 1 in binary slepian-wolf
    #p1g1_bsw = 0.85                     #Probability of 1 given 1 in binary LDPC
    #b_swbits_per_obit = 16              #Number of original sequence bits in each bit of the binary slepian-wolf code
    #p1g1_nbs = 0.7            #The probability of "coincidence" in nonbinary code
    #nb_swbits_per_obit = 64   #Effective number of original sequence bits per nonbinary letter.


    #Stats for the actual implementation of the codes
    #b_testparity = 965/1000.0       #Parities I need to send in binary code

    #nb_testparity = 540/1000.0      #Parities I need to send in non-binary code divided by total number of non-parity letters


    print "TOTAL SEQUENCE THEORETICAL ENTROPY:"
    ideal_entropy = theoretical(p1a,p1g1a,p1a,p1g1a)
    ideal_entropy2 = theoretical(p1b,p1g1b,p1b,p1g1b)
    print "\n\nFRAME OCCUPANCY THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is per alphabet bits of the original sequence
    b_entropy = theoretical_mat_nb(b_plet,b_mat,alph)/alph

    print "\nPercentage of total entropy recovered:",b_entropy/ideal_entropy2

    print "\n\nFRAME LOCATION THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is for a large amount of original sequence bits
    nb_entropy = theoretical_mat_nb(nb_plet,nb_mat,alph)/nb_bperf

    print "\nPercentage of total entropy recovered:",nb_entropy/ideal_entropy2

    print "\n\nTOTAL PERCENTAGE OF ENTROPY RETAINED:"
    print "Theory: %f"%((nb_entropy+b_entropy)/ideal_entropy)
    print "Theory2:%f"%((nb_entropy+b_entropy)/ideal_entropy2)
    
    print "\n\n"
    
    return (ideal_entropy,ideal_entropy2,b_entropy,nb_entropy)
