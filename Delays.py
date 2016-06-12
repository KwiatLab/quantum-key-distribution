'''
Created on Jun 8, 2016

@author: laurynas
'''
from numpy import concatenate,take,uint64,save, int64
import sys
import graphs
import ttag
from CalculateDelays import getDelays

def check_correlations(resolution, A_B_timetags, A_B_channels,channels1,channels2,d1,d2, coincidence_window_radius):
    print("- Applying Delays")
    
    for i in range(len(channels1)):
        print "-> ",channels1[i],A_B_timetags[A_B_channels==channels1[i]]," - ", d1[i]
        A_B_timetags[A_B_channels==channels1[i]]-=d1[i]
    print "--------"

    for i in range(len(channels2)):
        print "-> ",channels2[i],A_B_timetags[A_B_channels==channels2[i]]," - ", d2[i]
        A_B_timetags[A_B_channels==channels2[i]]-=d2[i]
            
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)      
    print "SORTED TTAGS: ", A_B_timetags,A_B_channels
    buf_num = ttag.getfreebuffer()
    bufDelays = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    bufDelays.resolution = resolution
    bufDelays.channels = max(A_B_channels)+1
    bufDelays.addarray(A_B_channels,A_B_timetags.astype(uint64))
    length_in_bins = int(coincidence_window_radius/resolution)*2
    print "# of bins before plotting", length_in_bins
    print "time", (A_B_timetags[-1]-1)*resolution
    graphs.plotABCorrelations(bufDelays,channels1,channels2,pulsebin = resolution, time=(A_B_timetags[-1]-1)*resolution, bins = int(coincidence_window_radius/resolution)*2 )
    
def calculate_delays(aliceTtags,aliceChannels,bobTtags,bobChannels,
                    resolution= 156.25e-12,
                    coincidence_window_radius = 1.9e-7,
                    channel_array = [2,3,4,5]):

    A_B_timetags = concatenate([aliceTtags,bobTtags])
    A_B_channels = concatenate([aliceChannels,bobChannels])

    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)

    
    buf_num = ttag.getfreebuffer()
    bufN = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    bufN.resolution = resolution
    print max(A_B_channels)
    bufN.channels = max(A_B_channels)+1
    bufN.addarray(A_B_channels,A_B_timetags)
    
    
    channels1 = channels2 = channel_array
#    1.9e-7 biggest u can make and still get correlations this corresponds to 1458 bins in diameter of coincidence window
#    UPDATE: actaully you can take smaller fraction of the strings to determine delays but then you need to increase coincidence window
    (d1,d2) = getDelays(bufN,channels1,channels2,delaymax=coincidence_window_radius,time=(A_B_timetags[-1])-1*bufN.resolution)
    
#     d1,d2 = getDelays(bufN,channels1,channels2,delays1=d1,delays2=d2,delaymax=bufN.resolution*100,time=(A_B_timetags[-1]-1)*bufN.resolution)
#     print (d1/buf.resolution,d2/buf.resolution)
    print (d1/bufN.resolution,d2/bufN.resolution)

    d1 = (d1/bufN.resolution).astype(int64)
    d2 = (d2/bufN.resolution).astype(int64)
    print "_____>",d1,d2
    sys.stdout.flush()
    
    print "will now plotting corellations to check if it looks good."
    check_correlations(resolution, A_B_timetags.astype(int64), A_B_channels, channels1, channels2, d1, d2, coincidence_window_radius)
    print("Saving Alice and Bob delays to file.")
    save("./resultsLaurynas/Delays/aliceDelay.npy",d1)
    save("./resultsLaurynas/Delays/bobDelay.npy",d2)
    

