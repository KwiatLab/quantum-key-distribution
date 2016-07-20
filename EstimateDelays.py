'''
Created on Jul 20, 2016

@author: laurynas
'''


from Delays import calculate_delays
from System import load_data
from ttag_delays import getDelay
from numpy import concatenate,take,uint64,uint8,save,zeros
import graphs
import ttag


def check_correlations(aliceTtags,aliceChannels,bobTtags,bobChannels,resolution, A_B_timetags, A_B_channels,channels1,channels2,delays,coincidence_window_radius,matrix_before,delay_max,dic,l):

    for delay,ch1,ch2 in zip(delays,channels1,channels2):
        if delay < 0:
            A_B_timetags[A_B_channels == ch2] += (abs(delay)).astype(uint64)
        else:
            A_B_timetags[A_B_channels == ch1] += delay.astype(uint64)
               
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)      




    buf_num = ttag.getfreebuffer()
    bufDelays = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    bufDelays.resolution = resolution
    bufDelays.channels = max(A_B_channels)+1
    bufDelays.addarray(A_B_channels,A_B_timetags.astype(uint64))
     

    A_B_timetags1 = concatenate([aliceTtags,bobTtags])
    A_B_channels1 = concatenate([aliceChannels,bobChannels])
 
    indexes_of_order = A_B_timetags1.argsort(kind = "mergesort")
    A_B_channels1 = take(A_B_channels1,indexes_of_order)
    A_B_timetags1 = take(A_B_timetags1,indexes_of_order)
  
    
    buf_num = ttag.getfreebuffer()
    bufCorrec = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    bufCorrec.resolution = resolution
    bufCorrec.channels = max(A_B_channels1)+1
    bufCorrec.addarray(A_B_channels1,A_B_timetags1.astype(uint64))



    length_in_bins = int(delay_max/resolution)*2
    
    graphs.plotABCorrelations(bufDelays,channels1,channels2,pulsebin = resolution, time=(A_B_timetags[-1]-1)*resolution, bins = length_in_bins)
    
def calculate_delays(aliceTtags,aliceChannels,bobTtags,bobChannels,
                    resolution= 78.125e-12,
                    coincidence_window_radius = 1000-12,
                    delay_max = 1e-7,
                    channels1, channels2):
    
    A_B_timetags = concatenate([aliceTtags,bobTtags])
    A_B_channels = concatenate([aliceChannels,bobChannels])

    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)

    buf_num = ttag.getfreebuffer()
    bufN = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e11))
    bufN.resolution = resolution
    bufN.channels = max(A_B_channels)+1
    bufN.addarray(A_B_channels,A_B_timetags)


    delays = zeros(len(channels1)+len(channels2))
    for i,j in zip(channels1, channels2):
        delays[i] = getDelay(bufN,i,j,delaymax=delay_max,time=(A_B_timetags[-1]-1)*bufN.resolution)

    
    print "will now plotting corellations to check if it looks good."
    check_correlations(aliceTtags,aliceChannels,bobTtags,bobChannels,resolution, A_B_timetags.astype(uint64), A_B_channels, channels1, channels2,delays/bufN.resolution, coincidence_window_radius, delay_max)
    print("Saving delays to file.")
    save("./resultsLaurynas/Delays/delays.npy",delays/bufN.resolution)
    

if __name__ == '__main__':


    alice_channels = [0,1,2,3]
    bob_channels =   [4,5,6,7]
    coincidence_window_radius = 1e-9
    (aliceTtags,aliceChannels) = load_data("alice",alice_channels,data_factor=10)
    (bobTtags,bobChannels) = load_data("bob",bob_channels,data_factor=10)

    indexes_of_order = aliceTtags.argsort(kind = "mergesort")
    aliceChannels = take(aliceChannels,indexes_of_order)
    aliceTtags = take(aliceTtags,indexes_of_order)
    
    indexes_of_order = bobTtags.argsort(kind = "mergesort")
    bobChannels = take(bobChannels,indexes_of_order)
    bobTtags = take(bobTtags,indexes_of_order)
    
    aliceTtags = aliceTtags[:len(aliceTtags)]
    aliceChannels = aliceChannels[:len(aliceChannels)]
    
    bobTtags = bobTtags[:len(bobTtags)]
    bobChannels = bobChannels[:len(bobChannels)]


    calculate_delays(aliceTtags.astype(uint64), aliceChannels.astype(uint8),
                     bobTtags.astype(uint64), bobChannels.astype(uint8),
                     coincidence_window_radius = coincidence_window_radius,
                     channels1 = alice_channels,channels2 = bob_channels) 


