'''
Created on Jun 8, 2016

@author: laurynas
'''
from numpy import concatenate,take,uint64,save, int64,zeros,where
import sys
import graphs
import ttag
from ttag_delays import getDelay, getDelays

def check_correlations(resolution, A_B_timetags, A_B_channels,channels1,channels2,delays,coincidence_window_radius,matrix_before):
    print "TIMETAGS BEFORE DELAYS:", A_B_timetags,A_B_channels
    print("- Applying Delays")
    
    for delay,ch1,ch2 in zip(delays,channels1,channels2):
        if delay < 0:
            print "abs-->>",(abs(delay)).astype(int)
            A_B_timetags[A_B_channels == ch2] += (abs(delay)).astype(uint64)
        else:
            print "-->",delay
            A_B_timetags[A_B_channels == ch1] += delay.astype(uint64)
               
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)      
    print "SORTED TTAGS: ", A_B_timetags,A_B_channels
    buf_num = ttag.getfreebuffer() 
    buffer = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    buffer.resolution = 78.125e-12
    buffer.channels = max(A_B_channels)+1
    buffer.addarray(A_B_channels,A_B_timetags)
    
    
#     (d1,d2) = getDelays(buffer,channels1,channels2,delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*buffer.resolution)
#     A_B_timetags = A_B_timetags.astype(float)
#     for i in range(len(channels1)):
#         A_B_timetags[A_B_channels==channels1[i]]-=d1[i]
#     for i in range(len(channels2)):
#         A_B_timetags[A_B_channels==channels2[i]]-=d2[i]
#     A_B_timetags = A_B_timetags.astype(uint64)
    
    
    buf_num = ttag.getfreebuffer()
    bufDelays = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    bufDelays.resolution = resolution
    bufDelays.channels = max(A_B_channels)+1
    bufDelays.addarray(A_B_channels,A_B_timetags.astype(uint64))
    
    print "__WITH_DELAYS-->\n",(bufDelays.coincidences((A_B_timetags[-1]-1)*bufDelays.resolution, coincidence_window_radius))
    print "__DIFF___->>>>\n",matrix_before.astype(int64)-(bufDelays.coincidences((A_B_timetags[-1]-1)*bufDelays.resolution, coincidence_window_radius).astype(int64))
    
    length_in_bins = int(coincidence_window_radius/resolution)*2
    print "# of bins before plotting", length_in_bins
    print "time", (A_B_timetags[-1]-1)*resolution
    graphs.plotABCorrelations(bufDelays,channels1,channels2,pulsebin = resolution, time=(A_B_timetags[-1]-1)*resolution, bins = int(coincidence_window_radius/resolution)*2 )
    
def calculate_delays(aliceTtags,aliceChannels,bobTtags,bobChannels,
                    resolution= 78.125e-12,
                    coincidence_window_radius = 1e-7):
    
    channels1 = [0,1,2,3]
    channels2 = [4,5,6,7]
    
    A_B_timetags = concatenate([aliceTtags,bobTtags])
    A_B_channels = concatenate([aliceChannels,bobChannels])

    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)

    
    buf_num = ttag.getfreebuffer()
    bufN = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    bufN.resolution = resolution
    bufN.channels = max(A_B_channels)+1
    bufN.addarray(A_B_channels,A_B_timetags)
    
    print "__BEFORE DELAYS-->\n",(bufN.coincidences((A_B_timetags[-1]-1)*bufN.resolution, coincidence_window_radius))
    
#    1.9e-7 biggest u can make and still get correlations this corresponds to 1458 bins in diameter of coincidence window
#    UPDATE: actaully you can take smaller fraction of the strings to determine delays but then you need to increase coincidence window
    delays = zeros(len(channels1))
    k = 0
    for i,j in zip(channels1, channels2):
        delays[i] = getDelay(bufN,i,j,delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*bufN.resolution)
#         delays[i] = getDelay(bufN,i,j,delaymax=coincidence_window_radius,time=5.0)
        print delays[i]
        k+=1
    
    delays1=zeros((len(channels1),len(channels2)))
    delays2=zeros((len(channels1),len(channels2)))
    
    for j in range(len(delays1)):
        for i in range(len(delays2)):
            delays2[j][i] = getDelay(bufN,channels1[j],channels2[i], delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*bufN.resolution)
            
    #Next, set all of delays for channels1
    for j in range(len(delays2)):
        for i in range(len(delays1)):
            delays1[j][i] = getDelay(bufN,channels1[i],channels2[j], delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*bufN.resolution)
    
    print "DIRECT_____>",(delays/bufN.resolution)
    print "MATRIX for alice_____>",(delays1/bufN.resolution)
    print "MATRIX for bob_____>",(delays2/bufN.resolution)

    sys.stdout.flush()
    
    print "will now plotting corellations to check if it looks good."
    check_correlations(resolution, A_B_timetags.astype(uint64), A_B_channels, channels1, channels2,delays/bufN.resolution, coincidence_window_radius, (bufN.coincidences((A_B_timetags[-1]-1)*bufN.resolution, coincidence_window_radius)))
#     
    print("Saving delays to file.")
    save("./resultsLaurynas/Delays/delays.npy",delays/bufN.resolution)
    

