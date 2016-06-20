'''
Created on Jun 11, 2016

@author: laurynas
'''
import ttag
from numpy import uint8,concatenate,uint64,loadtxt,array,take,set_printoptions,zeros,where
from ttag_delays import getDelay,getDelays

if __name__ == '__main__':
    alice = loadtxt("./DataFiles/FakeDelays/AliceFake2.txt")
    bob = loadtxt("./DataFiles/FakeDelays/BobFake2_negative.txt")
    
    aliceTtags= array(alice[1],dtype = uint64)
    bobTtags=array(bob[1], dtype = uint64)

    aliceChannels=array(alice[0],dtype = uint8)
    bobChannels = array(bob[0], dtype = uint8)
    
    A_B_timetags = concatenate((aliceTtags,bobTtags))
    A_B_channels = concatenate((aliceChannels,bobChannels))
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)
    A_B_channels.astype(uint8)
    A_B_timetags.astype(uint64)

    
    buf_num = ttag.getfreebuffer() 
    buffer = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    buffer.resolution = 156.25e-12
    buffer.channels = max(A_B_channels)+1
    buffer.addarray(A_B_channels,A_B_timetags)
    resolution= 156.25e-12,
    channels1 = [0,1,2,3]
    channels2 = [4,5,6,7]

#     (d1,d2) = calculate_delays(aliceTtags, aliceChannels, bobTtags, bobChannels,resolution= 156.25e-12,coincidence_window_radius = 5e-9, channel_array=channels1) 
#     print d1,d2
    set_printoptions(edgeitems = 50)

    coincidence_window_radius =1e-4
    print "radius in bins: ", coincidence_window_radius/buffer.resolution, "time back in bins, ", int(((A_B_timetags[-1]-1)*buffer.resolution)/buffer.resolution)
    
    delays = zeros(len(channels1),dtype = uint64)
    k = 0
    for i,j in zip(channels1,channels2):
        delay= (getDelay(buffer,i,j,delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*buffer.resolution)/buffer.resolution).astype(uint64)
        delays[k] = delay
        k+=1
    
    print "__BEFORE DELAYS-->\n",(buffer.coincidences((A_B_timetags[-1]-1)*buffer.resolution, coincidence_window_radius))

    
    (d1,d2) = getDelays(buffer,channels1,channels2,delaymax=coincidence_window_radius,time=(A_B_timetags[-1]-1)*buffer.resolution)
    
    for i in range(len(channels1)):
        A_B_timetags[A_B_channels==channels1[i]]-=d1[i]
    for i in range(len(channels2)):
        A_B_timetags[A_B_channels==channels2[i]]-=d2[i]
    
    
#     for delay,ch1,ch2 in zip(delays,channels1,channels2):
#         if delay < 0:
#             A_B_timetags[A_B_channels == ch2] += (abs(delay)).astype(uint64)
#         else:
#             A_B_timetags[A_B_channels == ch1] += delay.astype(uint64)
            
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)
    
    
    buf_num = ttag.getfreebuffer()
    bufDelays = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    bufDelays.resolution = resolution
    bufDelays.channels = max(A_B_channels)+1
    bufDelays.addarray(A_B_channels,A_B_timetags.astype(uint64))
    
    print "__WITH_DELAYS-->\n",(bufDelays.coincidences((A_B_timetags[-1]-1)*bufDelays.resolution, coincidence_window_radius))
       
    