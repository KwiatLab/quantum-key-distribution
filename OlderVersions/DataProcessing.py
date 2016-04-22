'''
Created on Feb 17, 2016

@author: laurynas
'''

import ttag
import numpy
from numpy import *
import sys


def read_raw_file(alice_raw_filename, bob_raw_filename, resolution=None):
    alice_file_contents = loadtxt(alice_raw_filename)
    bob_file_contents = loadtxt(bob_raw_filename)
    
    indexes_of_order = alice_file_contents[1, :].argsort(kind="mergesort")
    alice_file_contents = take(alice_file_contents, indexes_of_order, 1)
    indexes_of_order = bob_file_contents[1, :].argsort(kind="mergesort")
    bob_file_contents = take(bob_file_contents, indexes_of_order, 1)
    
        
    alice__channels = alice_file_contents[0, :].astype(uint8);
    bob__channels = bob_file_contents[0, :].astype(uint8);

    if (resolution):
        alice_timetags = around((alice_file_contents[1, :] -
                                 alice_file_contents[1, 0]) / resolution).astype(uint64)
        bob_timetags = around((bob_file_contents[1, :] -
                                 bob_file_contents[1, 0]) / resolution).astype(uint64)
    else:
        alice_timetags = (alice_file_contents[1, :] -
                                alice_file_contents[1, 0]).astype(uint64)
        bob_timetags = (bob_file_contents[1, :] -
                                bob_file_contents[1, 0]).astype(uint64)
    print("alice channels")
    print(alice__channels)
    print("bob channels")
    print(bob__channels)
    print("alice ttags")
    print(alice_timetags)
    print("bob ttags")
    print(bob_timetags)
    
    return ([alice__channels,alice_timetags,bob__channels,bob_timetags])
def read_processed_file(alice_processed_filename, bob_processed_filename):

    alice_raw = load("./results/" + alice_processed_filename + ".npy")
    bob_raw = load("./results/" + bob_processed_filename + ".npy")
    return (alice_raw, bob_raw)

def saveprep(fname,alice,bob,alice_pol,bob_pol):
    print("Saving Data Arrays")
    sys.stdout.flush()
    save("./resultsLaurynas/alice_"+fname+".npy",alice)
    save("./resultsLaurynas/bob_"+fname+".npy",bob)
    save("./resultsLaurynas/alice_pol_"+fname+".npy",alice_pol)
    save("./resultsLaurynas/bob_pol_"+fname+".npy",bob_pol)
    
def binData(buf,binsize):
    channels,timetags = buf[:]
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

if (__name__ == '__main__'):
   # alice_raw_filename = raw_input("Insert filename with Alice counts: ")
    # bob_raw_filename = raw_input("Insert filename with Bob counts: ")
    numpy.set_printoptions(edgeitems = 100)

    alice_raw_filename = "power238_1s_0.csv"
    bob_raw_filename = "power238_1s_1.csv"
    
    resolution = 5e-11
    alice_bob_ch_ttag =read_raw_file(alice_raw_filename, bob_raw_filename,resolution)
    saveprep("06-03-2014",alice_bob_ch_ttag[1],alice_bob_ch_ttag[3],alice_bob_ch_ttag[0],alice_bob_ch_ttag[2])
    
    buf_num = ttag.getfreebuffer()
#     buf = ttag.TTBuffer(buf_num,create=True,datapoints = len(alice_bob_ch_ttag[0]))
#     buf.resolution = resolution
#     buf.channels = max(alice_bob_ch_ttag[0])+1
#     buf.addarray(alice_bob_ch_ttag[0],alice_bob_ch_ttag[1])
#     buf.addarray(alice_bob_ch_ttag[2],alice_bob_ch_ttag[3])
#     print(len(alice_bob_ch_ttag[0]),len(alice_bob_ch_ttag[1]),len(alice_bob_ch_ttag[2]),len(alice_bob_ch_ttag[3])
#     print (buf.coincidences(2.0, 1e-5))
    print("Combining arrays")
     
    A_B_timetags = concatenate((alice_bob_ch_ttag[1],alice_bob_ch_ttag[3]))
    A_B_channels = concatenate((alice_bob_ch_ttag[0],alice_bob_ch_ttag[2]))
    #print A_B_timetags
    #print A_B_channels
    
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)
    print A_B_channels
    print("alice and bob timetags")
    print A_B_timetags 
    buf_num = ttag.getfreebuffer()
    print("Opening Buffer",buf_num)
    buf = ttag.TTBuffer(buf_num,create=True,datapoints = len(A_B_channels))
    
    print("Setting Properties")
    buf.resolution = resolution
    buf.channels = max(A_B_channels)+1
    print("->Resolution:",buf.resolution)
    print("->Channels:",buf.channels)
    
    print("Adding Data to Buffer")
    buf.addarray(A_B_channels,A_B_timetags)
      
    
    #---------------------------copied from prep- prep()
#     aliceChannels=[0,1,2,3,4,5]
#     bobChannels=[0,1,2,3,4,5]
# 
#     #binsize: The time in seconds of a bin
#     binsize = 1/(32*120.0e6)
# 
#     #Create the buffer that will do stuff
#     #buf = ttag.TTBuffer(1)
# 
#     print("Binning...")
# 
#     #The time allocated to each bin:
#     c,t = binData(buf,binsize)
#     print("c")
#     print (c)
#     print("t")
#     print (t)
#     
#     bob,alice,bob_pol,alice_pol = extractAliceBob(c,t,aliceChannels,bobChannels)
#     
#     
#     print("Finding Intersect")
#     #Make sure Alice and Bob datasets coencide
#     aliceMask = logical_and(alice > bob[0],alice < bob[-1])
#     bobMask = logical_and(bob > alice[0],bob < alice[-1])
# 
#     bob = bob[bobMask]
#     bob_pol = bob_pol[bobMask]
#     alice = alice[aliceMask]
#     alice_pol = alice_pol[aliceMask]
# 
#     #Now, rezero
#     z = min(bob[0],alice[0])
# 
#     bob-=z
#     alice-=z
#     print("bob,          alice,      bob_pol,alice_pol")
#     print(bob,alice,bob_pol,alice_pol)
    #----------------------------------------
    #print("Buffer number: ",buf_num,"is ready.")
    
    
    print (buf.coincidences(1.8, 1e-10))
    #alice_processed_filename = "alice"
    #bob_processed_filename = "bob_main_high"
    
    
#     #reading .csv and converting to numpy array
#     read_raw_file(alice_raw_filename, bob_raw_filename)
#    (alice_raw,bob_raw) = read_processed_file(alice_processed_filename,bob_processed_filename)
#     numpy.set_printoptions(edgeitems = 100)
#     alice_raw/=16
#    print(alice_raw)
#     print(bob_raw)
#     print(alice_raw[-1])
#     print(size(alice_raw))
#     print(size(bob_raw))
#     print(bob_raw)
    
    
    
    
