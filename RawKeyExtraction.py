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
from DataProcessing import *
from Delays import calculate_delays
from System import loadprep
# def loadprep(fname):
#     # print("Loading Alice and Bob arrays")
# 
#     sys.stdout.flush()
#     alice = load("./resultsLaurynas/resultsLaurynas/aliceTtags_"+fname+".npy")
#     bob = load("./resultsLaurynas/resultsLaurynas/bobTtags_"+fname+".npy")
#     alice_pol=load("./resultsLaurynas/resultsLaurynas/aliceChannels_"+fname+".npy")
#     bob_pol = load("./resultsLaurynas/resultsLaurynas/bobChannels_"+fname+".npy")
# 
#     return (alice,bob,alice_pol,bob_pol)

if (__name__ == '__main__'):
#     loadedData = loadprep("06032014_maxpower")
    
    
#     resolution = 1.5625e-10
#     resolution = 156.25e-12
#     print loadedData[0][:100]
#     numpy.set_printoptions(edgeitems = 20)
#==================================TRIAL FOR DELAYS==================================================
#     alice_raw_filename = "./DataFiles/ShorterFiles/06032014_maxpower_268_0_trimmed.csv"
#     bob_raw_filename = "./DataFiles/ShorterFiles/06032014_maxpower_268_1_trimmed.csv"
#     (aliceChannels,aliceTtags) = read_raw_file(alice_raw_filename, "alice", resolution)
#     
#     (bobChannels,bobTtags) = read_raw_file(bob_raw_filename, "bob", resolution)
# 
#     
#     aliceTtags = aliceTtags[:20]
#     aliceChannels = aliceChannels[:20]
# 
#     bobTtags = bobTtags[:20]
#     bobChannels = bobChannels[:20]
# 
#     
#     print "A: ",aliceChannels
#     print "B: ",bobChannels
#     
#     print "A: ",aliceTtags
#     print "B: ",bobTtags 
#     buf_num = ttag.getfreebuffer()
#     bufAlice = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e4))
#     bufAlice.resolution = resolution
#     bufAlice.channels = 6
#     bufAlice.addarray(aliceChannels,aliceTtags)
#     print bufAlice.singles((aliceTtags[-1])*bufAlice.resolution)
# 
#     buf_num = ttag.getfreebuffer()
#     bufBob = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e4))
#     bufBob.resolution = resolution
#     bufBob.channels = 6
#     bufBob.addarray(bobChannels,bobTtags)
#     #if subtract one from ttag to get rid of the error doesnt count first ttag so this's correct
#     print bufBob.singles((bobTtags[-1])*bufBob.resolution) 
#     
#     
#     A_B_timetags = concatenate([aliceTtags,bobTtags])
#     A_B_channels = concatenate([aliceChannels,bobChannels])
# 
#     indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
#     A_B_channels = take(A_B_channels,indexes_of_order)
#     A_B_timetags = take(A_B_timetags,indexes_of_order)
#     print "A&B: ",A_B_channels
#     print "A&B: ",A_B_timetags
#     
#     buf_num = ttag.getfreebuffer()
#     buf = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e4))
#     buf.resolution = resolution
#     buf.channels = max(A_B_channels)+1
#     buf.addarray(A_B_channels,A_B_timetags)
# 
#     
#     channels2=channels1 = [2,3,4,5]
#     coincidence_window_radius = 1.9e-7
#     (d1,d2) = getDelays(buf,channels1,channels2,delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*buf.resolution)
#     print (d1/buf.resolution,d2/buf.resolution)
#     
# #     graphs.plotABCorrelations(buf,channels1,channels2)
#     d1 = (d1/buf.resolution).astype(uint64)
# 
#     d2 = (d2/buf.resolution).astype(uint64)
#     print("- Applying Delays")
#     for i in range(len(channels1)):
#         A_B_timetags[A_B_channels==channels1[i]]-=d1[i]
#     for i in range(len(channels2)):
#         A_B_timetags[A_B_channels==channels2[i]]-=d2[i]
#         
#     indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
#     A_B_channels = take(A_B_channels,indexes_of_order)
#     A_B_timetags = take(A_B_timetags,indexes_of_order)    
#     print "Printing timetags after applying delays and sorting: ", A_B_timetags,A_B_channels
#     
#     buf_num = ttag.getfreebuffer()
#     bufDelays = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e4))
#     bufDelays.resolution = resolution
#     bufDelays.channels = max(A_B_channels)+1
#     bufDelays.addarray(A_B_channels,A_B_timetags)
#     
#     graphs.plotABCorrelations(bufDelays,channels1,channels2)
# =========================================================================================================

#     #adding Alice data to buffer 
#     buf_num = ttag.getfreebuffer()
#     # print("Opening Buffer",buf_num)
#     bufAlice = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
#     # print("Setting Properties")
#     bufAlice.resolution = resolution
#     bufAlice.channels = max(loadedData[2])+1
#     # print("->Resolution:",bufAlice.resolution)
#     # print("->Channels:",bufAlice.channels)
    (aliceTtags,aliceChannels) = loadprep("alice")
    (bobTtags,bobChannels) = loadprep("bob")
    
    alice_channels = [0,1,2,3]
    bob_channels =   [4,5,6,7]
    
#      1.9e-7 biggest u can make and still get correlations this corresponds to 1458 bins in diameter of coincidence window
#    UPDATE: actaully you can take smaller fraction of the strings to determine delays but then you need to increase coincidence window
    

#     # make them of equal size
    if (len(aliceTtags) > len(bobTtags)):
        aliceTtags    = aliceTtags[:len(bobTtags)]
        aliceChannels = aliceChannels[:len(bobChannels)]
    else:
        bobTtags    = bobTtags[:len(aliceTtags)]
        bobChannels = bobChannels[:len(aliceChannels)]

    indexes_of_order = aliceTtags.argsort(kind = "mergesort")
    aliceChannels = take(aliceChannels,indexes_of_order)
    aliceTtags = take(aliceTtags,indexes_of_order)
    
    indexes_of_order = bobTtags.argsort(kind = "mergesort")
#     
    bobChannels = take(bobChannels,indexes_of_order)
    bobTtags = take(bobTtags,indexes_of_order)
    calculate_delays(aliceTtags.astype(uint64), aliceChannels.astype(uint8), bobTtags.astype(uint64), bobChannels.astype(uint8)) 
    
#     #------------------------------------
#     # print("Alice ready. Adding Alice Data to Buffer")
#     bufAlice.addarray(aliceChannels,aliceTtags)
#     # print ('->>>bufAlice[:]', bufAlice[:])
# 
#     #creating Bob buffer ---------------------------------------------------------------------------
#     # print "Creating Bob Buffer"
#     buf_num = ttag.getfreebuffer()
# 
#     # print("Opening Buffer",buf_num)
#     
#     bufBob = ttag.TTBuffer(buf_num,create=True,datapoints = len(loadedData[1]))
#  
#     # print("Setting Properties")
#     bufBob.resolution = resolution
#     bufBob.channels = max(loadedData[3])+1
#     # print("->Resolution:",bufBob.resolution)
#     # print("->Channels:",bufBob.channels)
#     
#     #sorting Bob just in case ttags----------
#     # print ("Sorting Bob Tags")
   
#     #------------------------------------
#     #adding Bob data to buffer
#     # print("Adding Bob Data to Buffer")
#     '''
#     TO DO: figuree out why this line messes eveyrhing up. Probably getFreeBuffer doesnt work properly
#     '''
#     bufBob.addarray(bobChannels,bobTtags)
#     
#     # print ('->>>bufAlice[:]', bufBob[:])
#     # print "--------------------WILL BE CALCULATING STATISTICS-------------------------------------------------"
#     # printing statistics -------------------------------------------------------------------------------
# #     print aliceTtags[:100]
#     calculateStatistics(aliceTtags,bobTtags,aliceChannels,bobChannels)
#     # print "-------------------END OF STATISTICS---------------------------------------------------------------"
#     #adding buffer with ABdata---------------------------------------------------------------------------
#     # print("Combining ALICE and BOB and adding to new buffer")
#     # print aliceTtags.dtype
#     A_B_timetags = concatenate([aliceTtags,bobTtags])
#     A_B_channels = concatenate([aliceChannels,bobChannels])
#     # print A_B_timetags
#     # print A_B_channels
#     

# 
#     indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
#     A_B_channels = take(A_B_channels,indexes_of_order)
#     A_B_timetags = take(A_B_timetags,indexes_of_order)

    
#     buf_num = ttag.getfreebuffer()
#     buf = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
#     buf.resolution = resolution
#     buf.channels = max(A_B_channels)+1
#     buf.addarray(A_B_channels,A_B_timetags)
#     
    
#     channels2=channels1 = [2,3,4,5]
#    1.9e-7 biggest u can make and still get correlations this corresponds to 1458 bins in diameter of coincidence window
#    UPDATE: actaully you can take smaller fraction of the strings to determine delays but then you need to increase coincidence window
#     coincidence_window_radius = 1.9e-7
    
#     THIS IS ONLY TO SAVE DELAYS TO FILE FOR "System.py"
     
#     (d1,d2) = getDelays(buf,channels1,channels2,delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*buf.resolution)
#     print (d1/buf.resolution,d2/buf.resolution)
#      
# #     graphs.plotABCorrelations(buf,channels1,channels2)
#     d1 = (d1/buf.resolution).astype(uint64)
#     d2 = (d2/buf.resolution).astype(uint64)
#     print("- Applying Delays")
#     for i in range(len(channels1)):
#         A_B_timetags[A_B_channels==channels1[i]]-=d1[i]
#         aliceTtags[aliceChannels== channels1[i]]-=d1[i]
# 
#     for i in range(len(channels2)):
#         A_B_timetags[A_B_channels==channels2[i]]-=d2[i]
#         aliceTtags[bobChannels == channels2[i]]-=d2[i]
#         
#     indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
#     A_B_channels = take(A_B_channels,indexes_of_order)
#     A_B_timetags = take(A_B_timetags,indexes_of_order)      
#     
#     
#     indexes_of_order = aliceTtags.argsort(kind = "mergesort")
#     aliceChannels = take(aliceChannels,indexes_of_order)
#     aliceTtags = take(aliceTtags,indexes_of_order)
# 
#     indexes_of_order = bobTtags.argsort(kind = "mergesort")
#     bobChannels = take(bobChannels,indexes_of_order)
#     bobTtags = take(bobTtags,indexes_of_order)  
#   
#     print "A: ",aliceTtags,"\n"
#     print "B: ",bobTtags,"\n"
#     print "A&B",A_B_timetags,"\n"
#     
#     indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
#     A_B_channels = take(A_B_channels,indexes_of_order)
#     A_B_timetags = take(A_B_timetags,indexes_of_order)    
# 
#     A_B_D_timetags = concatenate([aliceTtags,bobTtags])
#     A_B_D_channels = concatenate([aliceChannels,bobChannels])
#     indexes_of_order = A_B_D_timetags.argsort(kind = "mergesort")
#     A_B_D_channels = take(A_B_D_channels,indexes_of_order)
#     A_B_D_timetags = take(A_B_D_timetags,indexes_of_order)
# 
#     
#     buf_num = ttag.getfreebuffer()
#     bufAD = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
#     bufAD.resolution = resolution
#     bufAD.channels = max(A_B_D_channels)+1
#     bufAD.addarray(A_B_D_channels,A_B_D_timetags)
# 
#     buf_num = ttag.getfreebuffer()
#     bufDelays = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
#     bufDelays.resolution = resolution
#     bufDelays.channels = max(A_B_channels)+1
#     bufDelays.addarray(A_B_channels,A_B_timetags)
#      
#     graphs.plotABCorrelations(bufAD,channels1,channels2)


#     # Coincidences ----------------------------------------------------------------------------
#     # print("\nCoincidence MatrixA:")
#     # print (bufAlice.coincidences(0.1, 1e-5))
#     
#     # print("\nCoincidence MatrixB:")
#     # print (bufBob.coincidences(0.1, 1e-5))
#     # print("\nCoincidence MatrixAB:")
#     # print (buf.coincidences(0.1, 1e-5))
#     
#     #-----------------------------Delays part -------------------------------------------------------------
#     # print ("Calculating delays")
#      
#     channels1=[2,3,4,5]
#     channels2=[2,3,4,5]
#     coincidence_window_radius = 6e-8
#     binsize = 260.41e-12
# #     d = getPossibleInitialDelays(buf,0,6)
#     print "->>>>>",aliceChannels,bobChannels, aliceTtags*resolution
#     d1,d2 = getDelays(buf,channels1,channels2,delaymax=coincidence_window_radius,time=1.0)
#     print ("d1->>>>",d1)
#     print ("d2->>>>",d2)
#     # print("Second Round of Delay finding")
# #     d1,d2 = getDelays(buf,channels1,channels2,delays1=d1,delays2=d2,delaymax=buf.resolution*100)
# #     print ("d1->>>>",d1)
# #     print ("d2->>>>",d2)
#     print("Preparing Correlation Plot")
#     bins = int(coincidence_window_radius/bufAlice.resolution)*2
#     print "bins"
# #   OLD  graphs.plotABCorrelations(buf,channels1,channels2,d1,d2)
#     graphs.plotABCorrelations(buf,channels1,channels2,d1,d2,bins,pulsebin=binsize)
# #     buf.__del__()
#     user=input("Looks good? (y/n):")
#     buf.__del__()
#     if (user=="y"):
#         print("Creating Syncd Data...")
#         channels,timetags = buf[:]
#         print "Printing timetags before: ", timetags
#     print("- Applying Delays")
#     for i in range(len(channels1)):
#         timetags[channels==channels1[i]]-=d1[i]
#     for i in range(len(channels2)):
#         timetags[channels==channels2[i]]-=d2[i]
#     print "Printing timetags after applying delays: ", timetags
#     
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
