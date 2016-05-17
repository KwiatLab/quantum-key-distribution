'''
Created on Mar 1, 2016

@author: laurynas
'''
""" The UIUC/NCSA license:

Copyright (c) 2014 Kwiat Quantum Information Group
All rights reserved.

Developed by:    Kwiat Quantum Information Group
                University of Illinois, Urbana-Champaign (UIUC)
                http://research.physics.illinois.edu/QI/Photonics/

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal with the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers
in the documentation and/or other materials provided with the distribution.

Neither the names of Kwiat Quantum Information Group, UIUC, nor the names of its contributors may be used to endorse
or promote products derived from this Software without specific prior written permission.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
"""

"""

"""

from numpy import *
import numpy
from matplotlib import *
from multiprocessing import Pool
from scipy.optimize import curve_fit
import ttag
import sys
from scipy.weave import inline
#import graphs
from SW_prep import *
from SlepianWolf import *
from entropy_calculator import *
from itertools import *

def probLetter(l,alph):
    p=zeros(alph)
    for i in xrange(alph):
        p[i]=sum(l==i)
    p/=len(l)
    return p

#Assumes that the person's data is ordered

def create_binary_string_from_laser_pulses(timetags, binsize = 263.15, relative_unit = 78.125):

    # change to int if possible (binsize / relative unit because laser frequency is 3.8HGz)
    BINSIZE = around(binsize/relative_unit)
    number_of_timetags = len(timetags)
    print "Number of ttags", number_of_timetags
    bin_string = zeros(number_of_timetags,dtype=uint64)


    for i in range(number_of_timetags):
        bin_number = around(timetags[i] / BINSIZE)
        bin_string[i]+=bin_number

    '''
    NOTE: If performance is really bad try fixing code below for implementation in low-level language
    '''
# 
#     code = """
#         long long z = 0;
#         long long bin_number = 0;
#         for (;z<number_of_timetags;z++) {
#             bin_number = round(timetags[z] / BINSIZE);
#             bin_string[bin_number] +=1;
#         }
#     """
#     inline(code,["BINSIZE","bin_string","timetags","number_of_timetags"],headers=["<math.h>"])

    return bin_string

def calculate_frame_occupancy(binary_string, frame_size):

    number_of_frames = around(binary_string[-1]/frame_size)+1
    print "Total number of frames: ", number_of_frames
    frame_occupancy = zeros(number_of_frames,dtype=uint16)
    
    for event_index in binary_string:
        frame_occupancy[event_index/frame_size] +=1    
        
    return frame_occupancy

'''
TO DO: Test manually first!!!! NEW: Passes simple tests
'''
def calculate_frame_locations(binary_string, frame_occupancies, frame_size):
    number_of_frames = len(frame_occupancies)
    frame_locations = zeros(number_of_frames,dtype=uint32)
    iterator = binary_string.__iter__()
    i=-1
    for element in iterator:
        i+=1
        # print "------------new iteration------------------"
        map_value = 0
        frame_number = element/frame_size
        # print "Frame number is: ", frame_number 
        position_in_frame = element%frame_size
        # print "Position in frame is: ", position_in_frame
        binary_position = frame_size - position_in_frame - 1
        # print "Binary position: ", binary_position
        map_value +=2**binary_position
        # print "Map value: ", map_value
        # to iterate through remaining elements in the frame
        for j in range (position_in_frame+1, frame_size):
            # print "position of element to be checked: ",(frame_number*frame_size+j)
            if (frame_number*frame_size+j) in binary_string  :
                # print "more elements to find"
                # print "degree of remaining",frame_size-j-1
                map_value +=2**(frame_size-j-1)
                iterator.next()
            else:
                continue
        frame_locations[frame_number] = map_value
        # print "Final map value: ", map_value
    return frame_locations

def createLDPCdata(timetags,polarizations,total_number_of_frames=None,frame_size=16):
    frame_numbers = timetags.copy()
    frame_numbers /= frame_size
    #Get the number of time bins
    if (total_number_of_frames == None): total_number_of_frames = frame_numbers[-1]+1

    #Create an array for the frame occupancy and location sequences
    frame_occupancy = zeros(total_number_of_frames,dtype=uint32)
    frame_location = zeros(total_number_of_frames,dtype=uint64)
    person_p = zeros(total_number_of_frames,dtype=uint8)
    frame_size = int(frame_size)

    print("Frame numbers: ", frame_numbers)
    number_of_timetags = int(argmax(frame_numbers>total_number_of_frames))
    if (number_of_timetags==0 and frame_numbers[0] > total_number_of_frames):
        print "ERROR: TOT is smaller than minvalue!"
        return None
    if (number_of_timetags==0):
        number_of_timetags = int(len(timetags))
    
    print("number of timetags: ",number_of_timetags)
    maxtag = (timetags[number_of_timetags-1])
    print("total number of frames: ",total_number_of_frames)

    code = """
        long long z = 0;
        for (;z<number_of_timetags;z++) {
            frame_occupancy[frame_numbers[z]] +=1;
            frame_location[frame_numbers[z]] += (timetags[z] % frame_size)*(pow((double)frame_size,(double)(frame_occupancy[frame_numbers[z]]-1)));
            person_p[frame_numbers[z]] = polarizations[z];
        }
    """
    inline(code,["frame_occupancy","frame_location","person_p","frame_numbers","polarizations","timetags","number_of_timetags","frame_size"],headers=["<math.h>"])

    # frame location counts how many occurances there are in the frame
    """
    for z in xrange(len(timetags)):
        frame_occupancy[frame_numbers[z]] += 1
        frame_location[frame_numbers[z]] = (timetags[z] % frame_size) * frame_occupancy[frame_numbers[z]]*frame_size
        person_p[frame_numbers[z]] = polarizations[z]
    """
    return (frame_occupancy,frame_location,person_p,maxtag)


#Assumes the timetags's data is ordered
def createBdata(person,tot=None):
    
    sw_p=person.copy()
    
    if (tot==None): tot=sw_p[-1]+1
    
    person_sw=zeros(tot,dtype=bool)
    
    doublenum=0
    for z in xrange(len(person)):
        #There is no value here, so save the value
        person_sw[sw_p[z]]=True
    
    return person_sw

def bincoincidences(A,B,pm=26):
    A=A.copy()
    B=B.copy()
    
    #Subtract the smallest time bin
    if (A[0]<=B[0]):
        A[:]-=A[0]-pm
        B[:]-=A[0]-pm
    else:
        A[:]-=B[0]-pm
        B[:]-=B[0]-pm
    
    #Split into time bins
    A/=pm*2
    B/=pm*2
    
    #Find coincidences
    coincidence = 0
    
    return len(intersect1d(A,B))

def binData(bufAlice,binsize):
    channels,timetags = bufAlice[:]
    timetags +=binsize/2-timetags[0]
    timebins = around(timetags/binsize).astype(uint64)
    timebins -= timebins[0] #zero the time bins
    return (channels,timebins)

def extractAliceBob(c,t,aliceChannels,bobChannels):
    bobmask = (c==bobChannels[0])
    for i in range(1,len(bobChannels)):
        bobmask = logical_or(bobmask,c==bobChannels[i])
    alicemask = (c == aliceChannels[0])
    for i in range(1,len(aliceChannels)):
        alicemask = logical_or(alicemask,c==aliceChannels[i])

    bob = t[bobmask]
    alice = t[alicemask]

    bob_pol = c[bobmask]
    alice_pol  = c[alicemask]

    #Reset the polarization detectors
    for i in range(len(aliceChannels)):
        alice_pol[alice_pol==aliceChannels[i]]=i
    for i in range(len(bobChannels)):
        bob_pol[bob_pol==bobChannels[i]]=i

    return (bob,alice,bob_pol,alice_pol)


def prep():
    #Alice and Bob Channels
    aliceChannels=[0,1,2,3]
    bobChannels=[4,5,6,7]

    #binsize: The time in seconds of a bin
    binsize = 1/(32*120.0e6)

    #Create the buffer that will do stuff
    bufAlice = ttag.TTBuffer(1)

    print("Binning...")

    #The time allocated to each bin:
    c,t = binData(bufAlice,binsize)

    #Create data sequences for Alice and Bob
    bob,alice,bob_pol,alice_pol = extractAliceBob(c,t,aliceChannels,bobChannels)

    print("Finding Intersect")
    #Make sure Alice and Bob datasets coencide
    aliceMask = logical_and(alice > bob[0],alice < bob[-1])
    bobMask = logical_and(bob > alice[0],bob < alice[-1])

    bob = bob[bobMask]
    bob_pol = bob_pol[bobMask]
    alice = alice[aliceMask]
    alice_pol = alice_pol[aliceMask]

    #Now, rezero
    z = min(bob[0],alice[0])

    bob-=z
    alice-=z

    return (alice,bob,alice_pol,bob_pol)


def saveprep(fname,alice,bob,alice_pol,bob_pol):
    print("Saving Data Arrays")
    sys.stdout.flush()
    save("./results/results/alice_"+fname+".npy",alice)
    save("./results/results/bob_"+fname+".npy",bob)
    save("./results/results/alice_pol_"+fname+".npy",alice_pol)
    save("./results/results/bob_pol_"+fname+".npy",bob_pol)

def loadprep(fname):
    print("Loading Alice and Bob arrays")

    sys.stdout.flush()
    alice = load("./resultsLaurynas/aliceTtags_"+fname+".npy")
    bob = load("./resultsLaurynas/bobTtags_"+fname+".npy")
    alice_pol=load("./resultsLaurynas/aliceChannels_"+fname+".npy")
    bob_pol = load("./resultsLaurynas/bobChannels_"+fname+".npy")

    return (alice,bob,alice_pol,bob_pol)

def resequence1(num,alph,times,seq=None,j=0):
    if (seq==None):
        seq = zeros(alph,dtype=bool)
    for i in range(times):
        seq[j+num%alph] = True
        num=num/alph
    return seq

def resequence(num,tries,alphabet):
    seq = zeros(alphabet*len(num),dtype=bool)
    for i in xrange(len(num)):
        resequence1(num[i],alphabet,tries[i],seq,alphabet*i)
    return seq


def calculateStatistics(alice,bob,alice_pol,bob_pol):
    #saveprep("main_high",*prep())
    numpy.set_printoptions(edgeitems = 100) 
    
    
    # print("alice channels")
    # print(alice_pol)
    print("alice ttags")
   # print(alice)
    # print ("bob channels")
    # print(bob_pol)
    print("bob ttags")

    binary_string_laser = create_binary_string_from_laser_pulses(alice)

    #Create the LDPC arrays (range 1-13)
    for frame_size in 2**array(range(1,2)):
        print("DOING ALPHABET",frame_size)
        # totlen = 8000000#min(max(int(alice[-1]/16),int(bob[-1]/16)),500000000)
        print("Extracting Data")
        frame_occupancies = calculate_frame_occupancy(binary_string_laser, frame_size)
        
#         can test frame_location algorithm with following line
#         print calculate_frame_locations(array([0,1,2,3,4,5,6,7,8,9,10,11,12]), array([4,4,4,1]), 4)

        # (bob_fo,bob_fl,bob_p,maxtag_b) = createLDPCdata(bob,bob_pol,total_number_of_frames = totlen,frame_size=frame_size)
        # print("frame occupancy NEW")
        # print(alice_fo)
        # print("frame location sequnces NEW")
        # print(alice_fl)
        # print "->Extracting ideal FRAME_OCCUPANCY (1-1) and POLARIZATION arrays"
        # sys.stdout.flush()
        
        # #2-1,2-2,etc:
        # bigmask = logical_or(alice_fo>1,bob_fo>1)
        # bigalice = alice_fl[bigmask]
        # bigbob = bob_fl[bigmask]
        # bigalice_t = alice_fo[bigmask]
        # bigbob_t = bob_fo[bigmask]

        # multibob = resequence(bigbob,bigbob_t,frame_size)
        # multialice = resequence(bigalice,bigalice_t,frame_size)

        # #1-1,etc

        # nb_mask = logical_and(alice_fo==1,bob_fo==1)
    
        # alice_nbf = alice_fl[nb_mask]
        # bob_nbf = bob_fl[nb_mask]
    
        # alice_pf = alice_p[nb_mask]
        # bob_pf = bob_p[nb_mask]
        # print "DATA saved for LDPC procedure"
        # savetxt("FRAME_128_DATA.csv",(alice_nbf,bob_nbf),fmt="%i")

        # #graphs.polarizationPlot(alice_pf,bob_pf)

        # print(len(alice_nbf))
        # print(alice_nbf,bob_nbf)

        # #Extract code statistics


        # print "->Code statistics"
        # sys.stdout.flush()
    
        # print "Frame Occupancy:"
        # print "\tLength:",len(bob_fo)
        # b_co = sum(bob_fo==alice_fo)
        # b_cor = float(b_co)/len(bob_fo)
        # print "\tCoincidence:",b_co,b_cor
        # print "\tError:",1-b_cor
    
        # print "Frame Location:"
        # print "\tLength:",len(bob_nbf)
        # nb_co = sum(alice_nbf==bob_nbf)
        # nb_cor = float(nb_co)/len(bob_nbf)
        # print "\tCoincidence:",nb_co,nb_cor
        # print "\tError:",1-nb_cor



        # #The probability of a 1 in the original-original sequence
        # p1 = float(len(alice))/alice[-1]
        # print "p1:",p1,float(len(bob))/bob[-1]
    
        # total_c = intersect1d(alice,bob)
        # p1g1 = float(len(total_c))/len(alice)
        # print "Coincidence rate (p1g1):",p1g1,float(len(total_c))/len(bob)

        # if (any(alice_fo >= frame_size) or any(bob_fo >= frame_size)):
        #     print "WARNING: Over the TOP!"
        #     alice_fo[alice_fo >= frame_size]=frame_size-1
        #     bob_fo[bob_fo >= frame_size]=frame_size-1
        # swtransmat = transitionMatrix_data2(alice_fo,bob_fo,frame_size)
        # swpl = probLetter(alice_fo,frame_size)
        # print "Letter Probabilities:"
        # print swpl
        # print "Transition Matrix (SW):"
        # print swtransmat

        # nbtransmat = transitionMatrix_data2(alice_nbf,bob_nbf,frame_size)
        # nbpl = probLetter(alice_nbf,frame_size)
        # print "Letter Probabilities:"
        # print nbpl
        # print "Transition Matrix (NB):"
        # print nbtransmat

        # #The total number of original bins is alice, and the final bits is the number of nonbinary left
        # nb_bperf= maxtag_a/len(bob_nbf)
        # print "Number of original bits per nonbinary:",nb_bperf


        # #2-x theory
        # multi_c = float(sum(logical_and(multialice,multibob)))
        # multi_p1b = float(sum(multibob))/float(len(multibob))
        # multi_p1g1b = multi_c/float(sum(multibob))
        # multi_p1a = float(sum(multialice))/float(len(multialice))
        # if(sum(multialice)!=0):
        #     multi_p1g1a = multi_c/float(sum(multialice))
        # else:
        #     multi_p1g1a = multi_c/float(1)
        # multi_bperf = maxtag_a/float(len(multibob))
        # print "MULTI"
        # print "Length:",len(multibob)
        # print "Ones:",sum(multibob),sum(multialice)
        # print "Coincidences",multi_c
        # print "p1",multi_p1b,multi_p1a
        # print "p1g1",multi_p1g1b,multi_p1g1a
        # print "Number of original bits per multi:",multi_bperf
        # entropy_left = theoretical(multi_p1a,multi_p1g1a,multi_p1b,multi_p1g1b)
        # multientropy = entropy_left/multi_bperf


        # (te,te2,be,nbe)=entropy_calculate2(p1,p1g1,p1,0.27,frame_size,swpl,swtransmat,nb_bperf,nbpl,nbtransmat)
        # f=open("results/THEORY_2014_6_3_high.csv","a")
        # f.write(str(frame_size)+" "+str(te)+" "+str(te2)+" "+str(be)+" "+str(nbe)+" "+str(multientropy)+"\n")
        # f.close()