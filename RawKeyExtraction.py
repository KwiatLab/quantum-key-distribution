'''
Created on Feb 28, 2016

@author: laurynas
'''
import ttag
import numpy
from numpy import *
import sys
from CalculateDelays import *
from scipy.optimize import curve_fit
import graphs
import matplotlib
from matplotlib import *
from itertools import product
import Statistics
from Statistics import calculateStatistics

def loadprep(fname):
    # print("Loading Alice and Bob arrays")

    sys.stdout.flush()
    alice = load("./resultsLaurynas/resultsLaurynas/aliceTtags_"+fname+".npy")
    bob = load("./resultsLaurynas/resultsLaurynas/bobTtags_"+fname+".npy")
    alice_pol=load("./resultsLaurynas/resultsLaurynas/aliceChannels_"+fname+".npy")
    bob_pol = load("./resultsLaurynas/resultsLaurynas/bobChannels_"+fname+".npy")

    return (alice,bob,alice_pol,bob_pol)

if (__name__ == '__main__'):
    loadedData = loadprep("06032014_maxpower") 
    resolution = 5e-11
#     print loadedData[0][:100]
    numpy.set_printoptions(edgeitems = 100)

    #adding Alice data to buffer 
    buf_num = ttag.getfreebuffer()
    # print("Opening Buffer",buf_num)
    bufAlice = ttag.TTBuffer(buf_num,create=True,datapoints = len(loadedData[0]))
    # print("Setting Properties")
    bufAlice.resolution = resolution
    bufAlice.channels = max(loadedData[2])+1
    # print("->Resolution:",bufAlice.resolution)
    # print("->Channels:",bufAlice.channels)
    
    aliceTtags = loadedData[0]
    aliceChannels =loadedData[2]

    bobTtags = loadedData[1]   
    bobChannels =loadedData[3]

    # make them of equal size
    if (len(aliceTtags) > len(bobTtags)):
        aliceTtags    = aliceTtags[:len(bobTtags)]
        aliceChannels = aliceChannels[:len(bobChannels)]
    else:
        bobTtags    = bobTtags[:len(aliceTtags)]
        bobChannels = bobChannels[:len(aliceChannels)]

    indexes_of_order = aliceTtags.argsort(kind = "mergesort")
    aliceChannels = take(aliceChannels,indexes_of_order)
    aliceTtags = take(aliceTtags,indexes_of_order)
    #------------------------------------
    # print("Alice ready. Adding Alice Data to Buffer")
    bufAlice.addarray(aliceChannels,aliceTtags)
    # print ('->>>bufAlice[:]', bufAlice[:])

    #creating Bob buffer ---------------------------------------------------------------------------
    # print "Creating Bob Buffer"
    buf_num = ttag.getfreebuffer()

    # print("Opening Buffer",buf_num)

    bufBob = ttag.TTBuffer(buf_num,create=True,datapoints = len(loadedData[1]))
 
    # print("Setting Properties")
    bufBob.resolution = resolution
    bufBob.channels = max(loadedData[3])+1
    # print("->Resolution:",bufBob.resolution)
    # print("->Channels:",bufBob.channels)
    
    #sorting Bob just in case ttags----------
    # print ("Sorting Bob Tags")
    indexes_of_order = bobTtags.argsort(kind = "mergesort")
    
    bobChannels = take(bobChannels,indexes_of_order)
    bobTtags = take(bobTtags,indexes_of_order)
    
    #------------------------------------
    #adding Bob data to buffer
    # print("Adding Bob Data to Buffer")
    '''
    TO DO: figuree out why this line messes eveyrhing up. Probably getFreeBuffer doesnt work properly
    '''
#     bufBob.addarray(bobChannels,bobTtags)
    
    # print ('->>>bufAlice[:]', bufBob[:])
    # print "--------------------WILL BE CALCULATING STATISTICS-------------------------------------------------"
    # printing statistics -------------------------------------------------------------------------------
#     print aliceTtags[:100]
#     calculateStatistics(aliceTtags,bobTtags,aliceChannels,bobChannels)
    # print "-------------------END OF STATISTICS---------------------------------------------------------------"
    #adding buffer with ABdata---------------------------------------------------------------------------
    # print("Combining ALICE and BOB and adding to new buffer")
    # print aliceTtags.dtype
    A_B_timetags = concatenate([aliceTtags,bobTtags])
    A_B_channels = concatenate([aliceChannels,bobChannels])
    # print A_B_timetags
    # print A_B_channels
    
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)
    # print A_B_channels

    buf_num = ttag.getfreebuffer()
    print("Opening Buffer",buf_num)
    buf = ttag.TTBuffer(buf_num,create=True,datapoints = len(A_B_channels))
    
    # print("Setting Properties")
    buf.resolution = resolution
    buf.channels = max(A_B_channels)+1
    # print("->Resolution:",buf.resolution)
    # print("->Channels:",buf.channels)
    
    # print("Adding Data to Buffer")
    buf.addarray(A_B_channels,A_B_timetags)
    # Coincidences ----------------------------------------------------------------------------
    # print("\nCoincidence MatrixA:")
    # print (bufAlice.coincidences(0.1, 1e-5))
    
    # print("\nCoincidence MatrixB:")
    # print (bufBob.coincidences(0.1, 1e-5))
    # print("\nCoincidence MatrixAB:")
    # print (buf.coincidences(0.1, 1e-5))
    
    #-----------------------------Delays part -------------------------------------------------------------
    # print ("Calculating delays")
#     
#     channels1=[0,1,2,3,4,5]
#     channels2=[0,1,2,3,4,5]
#      
#     d1,d2 = getDelays(buf,channels1,channels2,time=1.5)
#     print ("d1->>>>",d1)
#     print ("d2->>>>",d2)
#     # print("Second Round of Delay finding")
#     d1,d2 = getDelays(buf,channels1,channels2,delays1=d1,delays2=d2,delaymax=buf.resolution*100)
#     print ("d1->>>>",d1)
#     print ("d2->>>>",d2)
#     print("Preparing Correlation Plot")
#     graphs.plotABCorrelations(buf,channels1,channels2,d1,d2)
#     user=input("Looks good? (y/n):")
#     if (user=="y"):
#         print("Creating Syncd Data...")
#         channels,timetags = buf[:]
#         print "Printing timetags before: ", timetags
#     print("- Applying Delays")
#     for i in range(len(channels1)):
#         timetags[channels==channels1[i]]-=d1[i]
#     for i in range(len(channels2)):
#         timetags[channels==channels2[i]]-=d2[i]
# #  
#     print "Printing timetags after applying delays: ", timetags
    
    # print("- Extracting Alice and Bob")
#     allWanted = (channels==channels1[0])
#     for i in range(1,len(channels1)):
#         allWanted= logical_or(allWanted,channels==channels1[i])
#     for c in channels2:
#         allWanted = logical_or(allWanted,channels==c)
# 
#     channels = channels[allWanted]
#     timetags = timetags[allWanted]
# 
#     c1b = []
#     c2b = []
#     for c in range(len(channels1)):
#         c1b.append(channels==channels1[c])
#     for c in range(len(channels2)):
#         c2b.append(channels==channels2[c])
# 
#     for i in range(len(channels1)):
#         channels[c1b[i]]=i
#     for i in range(len(channels2)):
#         channels[c2b[i]]=i+len(channels1)
#  
#         
#     #
#     ##WTF: This code causes a segfault later! I don't even... I don't have time right now to fix it.
#     #
    # print("- Finding intersect of data from both time taggers")
#     #Find the first and last time tags of the two time taggers
#     #   and then take only the intersecting sets
#     c1I = c1b[0]
#     for i in range(1,len(c1b)):
#         c1I= logical_or(c1I,c1b[i])
#     c2I = c2b[0]
#     for i in range(1,len(c2b)):
#         c2I= logical_or(c2I,c2b[i])
#      
#     ttmin = logical_and(timetags > min(timetags[c1I]),timetags > min(timetags[c2I]))
#     ttmax = logical_and(timetags < max(timetags[c1I]),timetags < max(timetags[c2I]))
#     tttot = logical_and(ttmin,ttmax)
#     timetags = timetags[tttot]
#     channels = channels[tttot]
#     
#  
    # print("- Sorting")
#     #Sort again to make sure everything is fine
#     order = timetags.argsort()
#     timetags = take(timetags,order)
#     channels = take(channels,order)
#  
    # print(len(channels),len(timetags))
#  
#  
    # print("- Creating Buffer")
#     buf_num = ttag.getfreebuffer()
#  
    # print("- Opening Buffer",buf_num)
#     buf2 = ttag.TTBuffer(buf_num,create=True,datapoints = len(channels))
#  
    # print("- Setting Properties")
#     buf2.resolution = bufAlice.resolution
#     buf2.channels = max(channels)+1
    # print("- > Resolution:",buf2.resolution)
    # print("- > Channels:",buf2.channels)
#  
    # print("- Converting timetags to BIN format")
#     #First: Make the smallest tag 0 to avoid possible negatives
#     timetags-=timetags[0]
#     #Convert to bins
#     timetags = (around((timetags)/buf2.resolution)).astype(uint64)
#  
    # print(timetags,channels)
    # print("- Adding to Buffer")
#     buf2.addarray(channels,timetags)
#  
    # print("Buffer",buf_num,"Ready.")
