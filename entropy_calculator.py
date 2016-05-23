# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 2012
"""
#from pylab import *
import numpy
def log2(x):
    if x > 0:
        return numpy.log2(x)
    return 0

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


def theoretical_entropy_transition_matrix(p_letter,transmat,alphabet):
    #p1g1 is coincidence
    #p0g1 is failcoincidence
    # print "Transition Matrix:",transmat
    # print "Alphabet:",alphabet
    
    #Each letter of the alphabet has the same probability
    #p_letter = 1/float(alphabet)
    
    p_letA = numpy.dot(p_letter,transmat.transpose())
#     print"---------->>"     
#     print("PA:",p_letA)
#     print("PB:",p_letter)
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
def theoretical(probability_of_one_in_a,coincidence_rate_a,probability_of_one_in_b,coincidence_rate_b):
    
    #Probability calculations
    probability_of_zero_in_a=1-probability_of_one_in_a
    probability_of_zero_in_b=1-probability_of_one_in_b
    coincidence_rate_zero_a = 1-coincidence_rate_a
    coincidence_rate_zero_b = 1-coincidence_rate_b
    
    p1g0a=coincidence_rate_zero_a*probability_of_one_in_a/probability_of_zero_in_a
    p1g0b=coincidence_rate_zero_b*probability_of_one_in_b/probability_of_zero_in_b
    p0g0a=1-p1g0a
    p0g0b=1-p1g0b
    
    print "A: p0   %f\tp1   %f\np0g1 %f\tp1g1 %f\np0g0 %f\tp1g0 %f\n"%(probability_of_zero_in_a,probability_of_one_in_a,coincidence_rate_zero_a,coincidence_rate_a,p0g0a,p1g0a)
    print "B: p0   %f\tp1   %f\np0g1 %f\tp1g1 %f\np0g0 %f\tp1g0 %f\n"%(probability_of_zero_in_b,probability_of_one_in_b,coincidence_rate_zero_b,coincidence_rate_b,p0g0b,p1g0b)
    
    #Entropy calculations
    entropy_per_bit_a = -probability_of_one_in_a*log2(probability_of_one_in_a)-probability_of_zero_in_a*log2(probability_of_zero_in_a)
    entropy_per_bit_b = -probability_of_one_in_b*log2(probability_of_one_in_b)-probability_of_zero_in_b*log2(probability_of_zero_in_b)
    
    print "A: Entropy Per Bit:",entropy_per_bit_a
    print "B: Entropy Per Bit:",entropy_per_bit_b
    
    entropy_g1a=-coincidence_rate_a*log2(coincidence_rate_a)-coincidence_rate_zero_a*log2(coincidence_rate_zero_a)
    entropy_g0a=-p1g0a*log2(p1g0a)-p0g0a*log2(p0g0a)
    
    entropy_agb = entropy_g1a*probability_of_one_in_a+entropy_g0a*probability_of_zero_in_a
    
    entropy_left = entropy_per_bit_a-entropy_agb
    if (entropy_left<=0.0):
        entropy_left =0.0

    print "Minimum Entropy Sent:",entropy_agb
    print "Maximum Entropy Left:",entropy_left
    
    if (numpy.isnan(entropy_left)): return 0.0
    return entropy_left



def entropy_calculate(

    probability_of_one_in_a=0.05,        #The probability of a photon in a time bin (p1)
    coincidence_rate_a=0.4,              #The data's heralding efficiency (coincidence rate) (p1 given 1)
    probability_of_one_in_b=0.05,        #The probability of a photon in a time bin (p1)
    coincidence_rate_b=0.4,              #The data's heralding efficiency (coincidence rate) (p1 given 1)
    alphabet_size=16,                    #The alphabet (frame size used)
    b_probability_of_one=0.05,           #Probability of 1 in binary code
    b_coincidence_rate_a=0.85,           #Coincidence rate of binary slepian-wolf 1s
    b_probability_of_one_in_b=0.05,      #Probability of 1 in binary code
    b_coincidence_rate_b=0.85,           #Coincidence rate of binary slepian-wolf 1s
    coincidence_rate_non_binary=0.7,     #Probability of coincidence in nonbinary code
    nb_bperf = 64,                       #Effective number of original sequence bits per nonbinary letter
    b_test = 1.0,                        #Percentage sent in parities
    nb_test = 0.8,                       #Percentage sent in parities nonbinary
    transition_matrix_non_binary=None    #Transition matrix for nonbinary code
    ):

    print "\n\nTHEORY:\n\n"

    print "TOTAL SEQUENCE THEORETICAL ENTROPY:"
    ideal_entropy = theoretical(probability_of_one_in_a,coincidence_rate_a,probability_of_one_in_b,coincidence_rate_b)

    print "\n\nSLEPIAN-WOLF THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is per alphabet bits of the original sequence
    b_entropy = theoretical(b_probability_of_one,b_coincidence_rate_a,b_probability_of_one_in_b,b_coincidence_rate_b)/alphabet_size

    print "\nPercentage of total entropy recovered:",b_entropy/ideal_entropy

    print "\n\nNON-BINARY THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is for a large amount of original sequence bits
    nb_entropy = 0.0
    
    if (transition_matrix_non_binary!=None):
        nb_entropy = theoretical_entropy_transition_matrix(transition_matrix_non_binary,alphabet_size)/nb_bperf
    else:
        nb_entropy = theoretical_nb(coincidence_rate_non_binary,alphabet_size)/nb_bperf

    print "\nPercentage of total entropy recovered:",nb_entropy/ideal_entropy

    #These are calculations of entropy retained for the actual codes
    test_b_entropy = (1-b_test)/alphabet_size
    test_nb_entropy = log2(alphabet_size)*(1-nb_test)/nb_bperf

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
    probability_of_one_in_a=0.05,       #The probability of a photon in a time bin (p1)
    coincidence_rate_a=0.4,             #The data's heralding efficiency (coincidence rate) (p1 given 1)
    probability_of_one_in_b=0.05,       #The probability of a photon in a time bin (p1)
    coincidence_rate_b=0.4,             #The data's heralding efficiency (coincidence rate) (p1 given 1)
    alph=16,                            #The alphabet (frame size used)
    binary_letter_probabilites=None,    #Binary Probability of each letter
    b_mat=None,                         #Transition matrix for binary code
    nb_bperf = 64,                      #Effective number of original sequence bits per nonbinary letter
    nonbinary_letter_probabilities=None,#nonbinary probability of each letter
    transition_matrix_non_binary=None   #Transition matrix for nonbinary code
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
    ideal_entropy_alice = theoretical(probability_of_one_in_a,coincidence_rate_a,probability_of_one_in_a,coincidence_rate_a)
    ideal_entropy_bob = theoretical(probability_of_one_in_b,coincidence_rate_b,probability_of_one_in_b,coincidence_rate_b)
    print "\n\nFRAME OCCUPANCY THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is per alphabet bits of the original sequence
    print "BINARY Letter PROB--->>"
    print binary_letter_probabilites
    binary_entropy = theoretical_entropy_transition_matrix(binary_letter_probabilites,b_mat,alph)/alph

    print "\nPercentage of total entropy recovered:",binary_entropy/ideal_entropy_bob

    print "\n\nFRAME LOCATION THEORETICAL (APPROX) ENTROPY"
    #Note that the entropy retained is for a large amount of original sequence bits
    nonbinary_entropy = theoretical_entropy_transition_matrix(nonbinary_letter_probabilities,transition_matrix_non_binary,alph)/nb_bperf

    print "\nPercentage of total entropy recovered:",nonbinary_entropy/ideal_entropy_bob

    print "\n\nTOTAL PERCENTAGE OF ENTROPY RETAINED:"
#     print "SHOULD BE LESS THAN ONE AND ARE:"
#     print "\t %f"%(nonbinary_entropy/ideal_entropy_alice)
    print "\t %f"%(binary_entropy/ideal_entropy_alice) 
    print "Theory Alice: %f"%((nonbinary_entropy+binary_entropy)/ideal_entropy_alice)
    print "Theory Bob:%f"%((nonbinary_entropy+binary_entropy)/ideal_entropy_bob)
    
    print "\n\n"
    
    return (ideal_entropy_alice,ideal_entropy_bob,binary_entropy,nonbinary_entropy)
