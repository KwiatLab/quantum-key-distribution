
'''
Created on Jun 8, 2016
@author: laurynas
'''
from numpy import concatenate,take,uint64,save, int64,zeros,where,argwhere,intersect1d,load,array,append,\
    savetxt,in1d
import sys
import graphs
import ttag
from ttag_delays import getDelay, getDelays
from Statistics import create_binary_string_from_laser_pulses as laser

def remake_coincidence_matrix(coincidence_matrix):
    channels = len(coincidence_matrix[0][:])
    width = channels/2
    height = len(coincidence_matrix[:][0])/2
    matrix = zeros((height,width))
    
    for i in range(height):
        for j in range(width):
            matrix[i][j] = coincidence_matrix[i][channels/2+j]
    return matrix

def get_coinc(alice_timetags,alice_channels,bob_timetags,bob_channels,ch1,ch2,window_radius):
    coinc = 0
    alice_ttags = alice_timetags[alice_channels == ch1]
    bob_ttags = bob_timetags[bob_channels == ch2]
    
    for a,b in zip(alice_ttags,bob_ttags):
        if b >= a-window_radius and b<=a+window_radius:
            coinc+=1
    return coinc

def check_correlations(aliceTtags,aliceChannels,bobTtags,bobChannels,resolution, A_B_timetags, A_B_channels,channels1,channels2,delays,coincidence_window_radius,matrix_before,delay_max,dic,l,b1,b2):
#     print "TIMETAGS BEFORE DELAYS:", A_B_timetags,A_B_channels
#     print("- Applying Delays")
#     save("./Debugging/aliceChD.npy",A_B_channels[in1d(A_B_channels,channels1)])
#     save("./Debugging/aliceTtagsD.npy",A_B_timetags[in1d(A_B_channels,channels1)])
#     
#     save("./Debugging/bobChD.npy",A_B_channels[in1d(A_B_channels,channels2)])
#     save("./Debugging/bobTtagsD.npy",A_B_timetags[in1d(A_B_channels,channels2)])
#     
    for delay,ch1,ch2 in zip(delays,channels1,channels2):
        if delay < 0:
#             print "abs-->>",(abs(delay)).astype(int)
            A_B_timetags[A_B_channels == ch2] += (abs(delay)).astype(uint64)
        else:
#             print "-->",delay
            A_B_timetags[A_B_channels == ch1] += delay.astype(uint64)
               
    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)      
#     print "SORTED TTAGS: ", A_B_timetags,A_B_channels


    buf_num = ttag.getfreebuffer() 
    buffer = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e7))
    buffer.resolution = 260e-12
    buffer.channels = max(A_B_channels)+1
    buffer.addarray(A_B_channels,A_B_timetags)

#     a_ttags = array([])
#     a_channels = array([])
#     for ch in [0,1,2,3]:
#         a_ttags = append(a_ttags, A_B_timetags[A_B_channels == ch])
#         a_channels = append(a_channels, A_B_channels[A_B_channels == ch])
#                 
#     indexes_of_order = a_ttags.argsort(kind = "mergesort")
#     a_channels = take(a_channels,indexes_of_order)
#     a_ttags = take(a_ttags,indexes_of_order)
#     
#     
#     
# 
#     a_laser_delay_string = laser(a_ttags[a_channels == 3], resolution, coincidence_window_radius)
# 
# 
# 
# 
# 
#     b_ttags = array([])
#     b_channels = array([])
#     for ch in [4,5,6,7]:
#         b_ttags = append(b_ttags, A_B_timetags[A_B_channels == ch])
#         b_channels = append(b_channels, A_B_channels[A_B_channels == ch])
#                 
#     indexes_of_order = b_ttags.argsort(kind = "mergesort")
#     b_channels = take(b_channels,indexes_of_order)
#     b_ttags = take(b_ttags,indexes_of_order)
# 
# 
#     savetxt("./DarpaQKD/a_ttags3-7",a_ttags[a_channels == 3],fmt='%15d')
# #     print "with getcoinc",get_coinc(a_ttags, a_channels, b_ttags, b_channels, 3, 7, int(coincidence_window_radius/resolution))
#     savetxt("./DarpaQKD/intersection3-7.npy",intersect1d(a_laser_delay_string,b_ttags[b_channels == 7]),fmt='%10d')
#     print "Coincidences WITH delays",len(intersect1d(a_laser_delay_string,b_ttags[b_channels == 7]))





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
     
    
    
#     print bufDelays.singles((A_B_timetags[-1]-1)*bufDelays.resolution)
    with_delays = (bufDelays.coincidences((A_B_timetags[-1]-1)*bufDelays.resolution, coincidence_window_radius))
#     print "__WITH_DELAYS-->\n",with_delays
    remade = remake_coincidence_matrix(with_delays)
    
    
#     print "========================================================================================================="
#     print coincidence_window_radius
#     print "__REMADE-->>\n",remade
    remade1 = array([  [(remade[0][0] + remade[0][1]+remade[1][0]+remade[1][1])/sum(sum(remade)), (remade[0][2] + remade[0][3]+remade[1][2]+remade[1][3])/sum(sum(remade))],
                      [(remade[2][0] + remade[2][1]+remade[3][0]+remade[3][1])/sum(sum(remade)), (remade[2][2] + remade[2][3]+remade[3][2]+remade[3][3])/sum(sum(remade))]
                  ])
    print remade[2][2],remade[3][3]
    dic[l]=(remade[2][2]+remade[3][3])-(b1+b2)
    return dic
#     print array([  [(remade[0][0] + remade[0][1]+remade[1][0]+remade[1][1]), (remade[0][2] + remade[0][3]+remade[1][2]+remade[1][3])],
#                   [(remade[2][0] + remade[2][1]+remade[3][0]+remade[3][1]), (remade[2][2] + remade[2][3]+remade[3][2]+remade[3][3])]
#               ])
    print "========================================================================================================="
# 
#     print "__CONTRIBUTION-->>\n", remade1
    
    save("./Debugging/aliceChD.npy",A_B_channels[in1d(A_B_channels,channels1)])
    save("./Debugging/aliceTtagsD.npy",A_B_timetags[in1d(A_B_channels,channels1)])
     
    save("./Debugging/bobChD.npy",A_B_channels[in1d(A_B_channels,channels2)])
    save("./Debugging/bobTtagsD.npy",A_B_timetags[in1d(A_B_channels,channels2)])

    
    
#     print"----",len(bobTtags)
#     bobTtags = load("./DarpaQKD/bobCorrectedT.npy")
#     bobChannels = load("./DarpaQKD/bobCorrectedC.npy")
#     aliceTtags = load("./DarpaQKD/aliceCorrectedT.npy")
#     aliceChannels = load("./DarpaQKD/aliceCorrectedC.npy")


#     a_laser_string = laser(aliceTtags[aliceChannels==3], resolution, coincidence_window_radius)
#     print a_laser_string
#     print "Coincidences WITH CORRECTIONS",len(intersect1d(a_laser_string,bobTtags[bobChannels == 7]))/float(len(bobTtags[bobChannels == 7]))



#     print len(bobTtags)
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

    with_corrections = (bufCorrec.coincidences((A_B_timetags1[-1]-1)*bufCorrec.resolution, coincidence_window_radius))
#     print "__WITH_Corrections-->\n",with_corrections
#     print "__REMADE-->>\n",remake_coincidence_matrix(with_corrections)
 

    
#     a_laser_string = laser(aliceTtags[:len(aliceTtags)/88000], resolution, int(coincidence_window_radius/resolution))
#     print "Coincidences with delays",len(intersect1d(a_laser_string,bobTtags[:len(aliceTtags)/88000]))
#     print "with delays: number of 0 channels and number of channel 4",len(A_B_timetags[A_B_channels == 1]),len(A_B_timetags[A_B_channels == 4])
#     print "COINcidences between 0 and 4 after delays: ", get_coinc(aliceTtags, aliceChannels, bobTtags, bobChannels, 0, 4, int(coincidence_window_radius/resolution))

#     print "__DIFF___->>>>\n",matrix_before.astype(int64)-(bufDelays.coincidences((A_B_timetags[-1]-1)*bufDelays.resolution, coincidence_window_radius).astype(int64))
    
    length_in_bins = int(delay_max/resolution)*2
#     print "# of bins before plotting", length_in_bins
#     print "time", (A_B_timetags[-1]-1)*resolution
#     graphs.plotABCorrelations(bufDelays,channels1,channels2,pulsebin = resolution, time=(A_B_timetags[-1]-1)*resolution, bins = length_in_bins)
    
def calculate_delays(aliceTtags,aliceChannels,bobTtags,bobChannels, l,dic,
                    resolution= 260e-12,
                    coincidence_window_radius = 200-12,
                    delay_max = 1e-6):
    
    channels1 = [0,1,2,3]
    channels2 = [4,5,6,7]
    
    A_B_timetags = concatenate([aliceTtags,bobTtags])
    A_B_channels = concatenate([aliceChannels,bobChannels])

    indexes_of_order = A_B_timetags.argsort(kind = "mergesort")
    A_B_channels = take(A_B_channels,indexes_of_order)
    A_B_timetags = take(A_B_timetags,indexes_of_order)

#     print "TIME: ", A_B_timetags[-1]*resolution, " in seconds"     
    buf_num = ttag.getfreebuffer()
    bufN = ttag.TTBuffer(buf_num,create=True,datapoints = int(5e11))
    bufN.resolution = resolution
    bufN.channels = max(A_B_channels)+1
    bufN.addarray(A_B_channels,A_B_timetags)
#     print aliceChannels[-20:],aliceTtags[-20:]
#     print bobChannels[-20:],bobTtags[-20:]
#     print A_B_channels[-20:],A_B_timetags[-20:]
    coincidences_before = (bufN.coincidences((A_B_timetags[-1]-1)*bufN.resolution, coincidence_window_radius))
#     print "__BEFORE DELAYS-->\n",coincidences_before
    remade = remake_coincidence_matrix(coincidences_before)
#     print "__REMADE-->>\n",remade
    remade1 = array([  [(remade[0][0] + remade[0][1]+remade[1][0]+remade[1][1])/sum(sum(remade)), (remade[0][2] + remade[0][3]+remade[1][2]+remade[1][3])/sum(sum(remade))],
                      [(remade[2][0] + remade[2][1]+remade[3][0]+remade[3][1])/sum(sum(remade)), (remade[2][2] + remade[2][3]+remade[3][2]+remade[3][3])/sum(sum(remade))]
                  ])
#     
#     print array([  [(remade[0][0] + remade[0][1]+remade[1][0]+remade[1][1]), (remade[0][2] + remade[0][3]+remade[1][2]+remade[1][3])],
#                   [(remade[2][0] + remade[2][1]+remade[3][0]+remade[3][1]), (remade[2][2] + remade[2][3]+remade[3][2]+remade[3][3])]
#               ])
    print "",remade[2][2],remade[3][3]
#     print "__CONTRIBUTION-->>\n", remade1
#     
#     
#     
#     
#     a_laser_string = laser(aliceTtags[aliceChannels==3], resolution, coincidence_window_radius)
#     print "Coincidences without delays",len(intersect1d(a_laser_string,bobTtags[bobChannels == 7]))/float(len(bobTtags[bobChannels == 7]))





#     print "number of 0 channels and number of channel 4",len(A_B_timetags[A_B_channels == 0]),len(A_B_timetags[A_B_channels == 4])
#     print "Coincidences between 0 and 4 by Laurynas: ", get_coincidences(A_B_timetags, A_B_channels, 0, 4, coincidence_window_radius*2/bufN.resolution)
#     print "That's it"
    
#    1.9e-7 biggest u can make and still get correlations this corresponds to 1458 bins in diameter of coincidence window
#    UPDATE: actaully you can take smaller fraction of the strings to determine delays but then you need to increase coincidence window
    delays = zeros(7)
    k = 0
    for i,j in zip(channels1, channels2):
        delays[i] = getDelay(bufN,i,j,delaymax=delay_max,time=(A_B_timetags[-1]-1)*bufN.resolution)
#         delays[i] = getDelay(bufN,i,j,delaymax=coincidence_window_radius,time=5.0)
#         print delays[i]
        k+=1
    
    delays1=zeros((len(channels1),len(channels2)))
    delays2=zeros((len(channels1),len(channels2)))
#     
#     for j in range(len(delays1)):
#         for i in range(len(delays2)):
#             delays2[j][i] = getDelay(bufN,channels1[j],channels2[i], delaymax=delay_max,time=(A_B_timetags[-1]-1)*bufN.resolution)
#             
#     #Next, set all of delays for channels1
#     for j in range(len(delays2)):
#         for i in range(len(delays1)):
#             delays1[j][i] = getDelay(bufN,channels1[i],channels2[j], delaymax=delay_max,time=(A_B_timetags[-1]-1)*bufN.resolution)
    delays1 = delays1/bufN.resolution
    delays2 = delays2/bufN.resolution
#     print "DIRECT_____>",(delays/bufN.resolution)
#     
#     print "MATRIX for alice_____>",(delays1)
#     print "MATRIX for bob_____>",(delays2)
#     
#     for i in range(len(channels1)):
#         for j in range(len(channels2)-1):
#             print channels2[j],"-",channels2[j+1],": ",delays1[i][j]-delays1[i][j+1] 

#     sys.stdout.flush()
    
#     print "will now plotting corellations to check if it looks good."
#     return check_correlations(aliceTtags,aliceChannels,bobTtags,bobChannels,resolution, A_B_timetags.astype(uint64), A_B_channels, channels1, channels2,delays/bufN.resolution, coincidence_window_radius,coincidences_before, delay_max,dic,l,remade[2][2],remade[3][3])
#     
#     print("Saving delays to file.")
    save("./resultsLaurynas/Delays/delays.npy",delays/bufN.resolution)
