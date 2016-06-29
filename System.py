'''
Created on Jun 4, 2016

@author: laurynas
'''
import sys
import threading
import json
from os import system
import DataProcessing
from numpy import *
import ttag
from multiprocessing import *
from subprocess import Popen, list2cmdline
from Statistics import *
from ParityCheckMatrixGen import gallager_matrix
from SlepianWolf import encode
import timeit

def get_coinc(alice_timetags,alice_channels,bob_timetags,bob_channels,ch1,ch2,window_radius):
    coinc = 0
    alice_ttags = alice_timetags[alice_channels == ch1]
    bob_ttags = bob_timetags[bob_channels == ch2]
    
    for a,b in zip(alice_ttags,bob_ttags):
        if b >= a-window_radius and b<=a+window_radius:
            coinc+=1
    return coinc

def get_timetag_corrections(timetags,resolution,sync_period,coincidence_window_radius):
    indexes = {}
    timetags = timetags.astype(uint64)
    
    sync_block_size = int(sync_period/resolution)
    D_block_size = coincidence_window_radius*2+1


    for i in range(len(timetags)):
        ith = (timetags[i] % sync_block_size) % D_block_size
        indexes[int(timetags[i]/sync_block_size)] = ith
    
    return indexes

def bob_corrected(bob_ttags,sync_period,resolution,coincidence_window_radius):
    bob_ttag_dict = {}
    sync_block_size = int(sync_period/resolution)
    D_block_size = coincidence_window_radius*2+1
    
    for i in range(len(bob_ttags)):
        ith = (bob_ttags[i] % sync_block_size)
        bob_ttag_dict[int(bob_ttags[i]/sync_block_size)] = ith
    return bob_ttag_dict

def do_correction(bob_ttag_dict, bob_dict,alice_dict, coincidence_window_radius):
    number_of_bins_in_block = coincidence_window_radius*2+1
    corrected_string = array([])
    
    for a_key in alice_dict.keys():
        if bob_dict.has_key(a_key):
            distance_away = min(number_of_bins_in_block-bob_dict[a_key]+alice_dict[a_key], abs(alice_dict[a_key] - bob_dict[a_key]), number_of_bins_in_block-alice_dict[a_key]+bob_dict[a_key])
#             print "\t",number_of_bins_in_block-bob_dict[a_key]+alice_dict[a_key], abs(alice_dict[a_key] - bob_dict[a_key]), number_of_bins_in_block-alice_dict[a_key]+bob_dict[a_key]
            if distance_away <= coincidence_window_radius:
                if alice_dict[a_key]-coincidence_window_radius > 0:
                    left_boundary = alice_dict[a_key]-coincidence_window_radius
                else:
                    left_boundary = number_of_bins_in_block + (alice_dict[a_key]-coincidence_window_radius)
       
                if alice_dict[a_key] + coincidence_window_radius <= number_of_bins_in_block:
                    right_boundary = alice_dict[a_key] + coincidence_window_radius
                else:
                    right_boundary = coincidence_window_radius - (number_of_bins_in_block - alice_dict[a_key])
    
#                 print alice_dict[a_key],distance_away,left_boundary,right_boundary
    
                shift = 0
                if bob_dict[a_key] >= left_boundary:
                    shift = distance_away
#                 elif bob_dict[a_key] < left_boundary:
#                     bob_dict[a_key] = -1
#                     continue
                elif bob_dict[a_key] <= right_boundary:
                    shift = -distance_away
#                 elif bob_dict[a_key] > right_boundary:
# #                     bob_dict[a_key] = -1
#                     continue

                bob_dict[a_key] += shift
                bob_ttag_dict[a_key]  += shift
            else:
                bob_dict[a_key] = -1

    return bob_ttag_dict
     
def make_equal_size(alice_thread,bob_thread):
    if (len(alice_thread.ttags) > len(bob_thread.ttags)):
        alice_thread.ttags = alice_thread.ttags[:bob_thread.ttags]
        alice_thread.channels = alice_thread.channels[:len(bob_thread.channels)]
    else:
        bob_thread.ttags    = bob_thread.ttags[:len(alice_thread.ttags)]
        bob_thread.channels = bob_thread.channels[:len(alice_thread.channels)]
def loadprep(name,channelArray):

    sys.stdout.flush()
    all_ttags = load("./DarpaQKD/"+name+"Ttags.npy")
    all_channels = load("./DarpaQKD/"+name+"Channels.npy")
    all_ttags = all_ttags[:len(all_ttags)/10]
    all_channels = all_channels[:len(all_channels)/10]
    
    ttags = array([])
    channels = array([])
    for ch in channelArray:
        ttags = append(ttags, all_ttags[all_channels == ch])
        channels = append(channels, all_channels[all_channels == ch])
                
    indexes_of_order = ttags.argsort(kind = "mergesort")
    channels = take(channels,indexes_of_order)
    ttags = take(ttags,indexes_of_order)
    
    return (ttags,channels)

def load_save_raw_file(dir, alice_channels, bob_channels):
    data = loadtxt(dir)

    channels = data[:,0]
    timetags = data[:,1]
    print("Saving Data Arrays")
    sys.stdout.flush()
    
    save("./DarpaQKD/aliceChannels.npy",channels[in1d(channels, alice_channels)])
    save("./DarpaQKD/aliceTtags.npy",timetags[in1d(channels, alice_channels)])
    
    save("./DarpaQKD/bobChannels.npy",channels[in1d(channels, bob_channels)])
    save("./DarpaQKD/bobTtags.npy",timetags[in1d(channels, bob_channels)])
    
    
def LDPC_encode(party,column_weight = 3,number_parity_edges = 6):
    total_string_length = len(party.received_string)
    
    number_of_parity_check_eqns_gallager = int(total_string_length*column_weight/number_parity_edges)
    party.parity_matrix = gallager_matrix(number_of_parity_check_eqns_gallager, total_string_length, column_weight, number_parity_edges)
    
    party.syndromes=encode(party.parity_matrix,party.non_zero_positions[:len(party.received_string)],party.frame_size)
 
def LDPC_decode(party,decoder='bp-fft', iterations=70, frozenFor=5):
    party.sent_string = party.non_zero_positions[:len(party.received_string)]
    
    transition_matrix = transitionMatrix_data2(party.sent_string,party.received_string,party.frame_size)
    prior_probability_matrix = sequenceProbMatrix(party.sent_string,transition_matrix)
    belief_propagation_system = SW_LDPC(party.parity_matrix, party.syndromes, prior_probability_matrix, original=alice_thread.sent_string,decoder=decoder)
    
    belief_propagation_system.decode(iterations=iterations,frozenFor=frozenFor)
    
class PartyThread(threading.Thread):
    def __init__(self, resolution,name,channelArray, coincidence_window_radius,delay_max,sync_period):
        threading.Thread.__init__(self)
        self.running = True
        self.name = name
        self.resolution = resolution
        self.raw_file_name = None
        self.delays = None
        self.channelArray = channelArray
        self.coincidence_window_radius = coincidence_window_radius
        self.sync_period = sync_period
        self.other_party_correction = None
        self.ttags = None
        self.channels = None
        self.buffer = None
        self.event = threading.Event()
        self.frame_size = 2
        self.binary_laser_string = None
        self.frame_occupancies = None
        self.frame_locations = None
        self.non_zero_positions = None
        self.received_string = None
        self.error_rate = None
        self.race_flag = False
        self.parity_matrix = None
        self.syndromes = None
        self.sent_string = None
        self.delay_max = delay_max
        self.event.set()
    def do_clear(self):
        self.event.clear()
    def do_set(self):
        self.event.set()
    def run(self):
        while self.running:
#             print self.name.upper()+" : Reading.csv files and converting to .npy\n"
#             system("python ./DataProcessing.py "+self.raw_file_name+" "+self.name)
            print self.name.upper()+": Loading .npy data"
#             print self.ttags
            (self.ttags,self.channels) = loadprep(self.name, self.channelArray)
#             print self.name, self.channels,self.ttags
#           TODO: Add delays here
            print "Loading delays"
            self.delays = load("./resultsLaurynas/Delays/delays.npy")
            print self.delays
            
            self.ttags=self.ttags.astype(int64)
            print "Applying Delays"
            
            print "BEFORE DELAYS", self.ttags,self.channels
            for delay,ch in zip(self.delays,self.channelArray):
                if delay < 0 and self.name == "bob":
                    print "abs-->>",(abs(delay)).astype(uint64)
                    self.ttags[self.channels == ch] += (abs(delay)).astype(uint64)
                elif delay >= 0 and self.name == "alice":
                    print "-->",delay
                    self.ttags[self.channels == ch] += delay.astype(uint64)
             

            
            indexes_of_order = self.ttags.argsort(kind = "mergesort")
            self.channels = take(self.channels,indexes_of_order)
            self.ttags = take(self.ttags,indexes_of_order)
#             print"--->>>>>>>>>>>>>>>",self.ttags
            self.ttags = self.ttags.astype(uint64)
            self.channels = self.channels.astype(uint8)
            
            print self.name.upper() +" FINISHED with data. Will notify main.\n"        
            buf_num = ttag.getfreebuffer() 
            self.buffer = ttag.TTBuffer(buf_num,create=True,datapoints = len(self.ttags))
            self.buffer.resolution = self.resolution
            self.buffer.channels = max(self.channels)+1
#             print self.channels
         
#             indexes_of_order = self.ttags.argsort(kind = "quicksort")
#             self.channels = take(self.channels,indexes_of_order)
#             self.ttags = take(self.ttags,indexes_of_order)
            
            #-----------------TODO: DEBUG ADD ARRAY-------------------
            # print("Alice ready. Adding Alice Data to Buffer")
#             print self.name+"Will try to add array to buffer"
#             self.buffer.addarray(self.channels,self.ttags)
            print self.name.upper()+": Waiting for OPTIMAL FR SIZE" 
# ==========TYPICAL BLOCK TO WAIT FOR MAIN and after RESET SELF AGAIN
            self.do_clear()
            while main_event.is_set():
                pass
            self.race_flag = True
# ===================================================================   
#             print "BEFORE ", self.ttags
            if self.name == "alice":
                print "!!!!!!!!!!"        
                self.binary_string_laser = get_timetag_corrections(self.ttags[self.channels == 0], self.resolution, self.sync_period, int(self.coincidence_window_radius/self.resolution))
                print self.ttags
#   print "----",self.binary_string_laser
                bob_thread.other_party_correction = self.binary_string_laser
            else:
                self.binary_string_laser = get_timetag_corrections(self.ttags[self.channels == 4], self.resolution, self.sync_period, int(self.coincidence_window_radius/self.resolution))
                print self.ttags
#             self.binary_string_laser = self.ttags
# #             print "New Binary laser string", self.binary_string_laser
#             print self.name+" Calculating frame occupancies and locations:\n"
#             self.frame_occupancies = calculate_frame_occupancy(self.binary_string_laser,self.frame_size)
#             self.frame_locations = calculate_frame_locations_daniels_mapping(self.binary_string_laser, self.frame_occupancies, self.frame_size)
#             print "Locations",self.name,"--->\n",self.frame_locations
#             print self.name.upper()+": Notifying main that finished with laser strings. Will be waiting for"

#===========READY TO ANNOUNCE======================
            self.do_clear()
            print "Done with laser string"

            while main_event.is_set():
                pass
            print self.name.upper()+": I was released and will do error calc.\n"    
#             print self.received_string,self.non_zero_positions[:len(self.received_string)]
            self.sent_string = self.non_zero_positions[:len(self.received_string)]
            self.error_rate = 1-float(sum(self.received_string == self.sent_string))/len(self.received_string)
            print self.name.upper() + " error rate: ",self.error_rate,"\n"
                
            
            self.running = False
               
if __name__ == '__main__':

    
#     alice_raw_filename = "./DataFiles/ShorterFiles/06032014_maxpower_268_0_trimmed.csv"
#     bob_raw_filename = "./DataFiles/ShorterFiles/06032014_maxpower_268_1_trimmed.csv"
    raw_file_dir = "./DarpaQKD/Alice1_Bob1.csv"
    alice_channels = [0,1,2,3]
    bob_channels =   [4,5,6,7]
    
#     load_save_raw_file(raw_file_dir, alice_channels, bob_channels)
    
    
    set_printoptions(edgeitems = 20)
    resolution = 78.125e-12
    coincidence_window_radius = 1200e-12
    delay_max = 1e-5
    sync_period = 7.8125e-9
    
    alice_event = threading.Event()
    alice_event.set()
   
    bob_event = threading.Event()
    bob_event.set() 
    
    alice_thread = PartyThread(resolution, "alice",channelArray=alice_channels, coincidence_window_radius = coincidence_window_radius,delay_max = delay_max,sync_period=sync_period)
    bob_thread = PartyThread(resolution,"bob", channelArray=bob_channels, coincidence_window_radius = coincidence_window_radius,delay_max = delay_max,sync_period=sync_period)
    start = timeit.default_timer()  
   
    main_event = threading.Event()
    main_event.set()
    alice_thread.start()
    bob_thread.start()

    print "MAIN: will wait till AB finished loading data."
    while alice_thread.event.is_set() or bob_thread.event.is_set():
        pass
#     print "BEFORE EQUAL"
#     print alice_thread.ttags
#     print bob_thread.ttags
#     print "MAIN: Making datasets equal."
#     make_equal_size(alice_thread, bob_thread)

    print "MAIN: STATISTICS: "
    (alice,bob,alice_chan,bob_chan) = (alice_thread.ttags, bob_thread.ttags, alice_thread.channels, bob_thread.channels)
    
    
#     statistics = calculateStatistics(alice,bob,alice_chan,bob_chan, laser_jitter, resolution)
#     print statistics
    
#     max_shared_binary_entropy = max(statistics.values())
#     optimal_frame_size = int(list(statistics.keys())[list(statistics.values()).index(max_shared_binary_entropy)])

    optimal_frame_size = 2048
    alice_thread.frame_size = optimal_frame_size
    bob_thread.frame_size = optimal_frame_size

    print "MAIN: Optimal size calculated and set for both threads, release THEM!"
#     alice_binary_laser_string = load("./DarpaQKD/aliceLaserString.npy")
#     bob_binary_laser_string = load("./DarpaQKD/bobLaserString.npy")    
#     print "-->>",alice_binary_laser_string[alice_chan == 0]
#     print "-->>",bob_binary_laser_string[bob_chan == 4]
#     print "coincide all",len(intersect1d(alice_binary_laser_string[alice_chan == 0], bob_binary_laser_string[bob_chan == 4])) ,float(len(intersect1d(alice_binary_laser_string[alice_chan == 0], bob_binary_laser_string[bob_chan == 4])))/len(alice_binary_laser_string[alice_chan == 0])
#     print ":-: ",intersect1d(alice_binary_laser_string[alice_chan == 0], bob_binary_laser_string[bob_chan == 4])
    main_event.clear()
    alice_thread.do_set()
    bob_thread.do_set()

    while not(alice_thread.race_flag and bob_thread.race_flag):
        pass
    main_event.set()

#   ================TYPICAL BLOCK TO WAIT FOR BOTH AND THEN RELEASE AND RESET FLAG AGAIN  
# ========================================================================================
#     print "The maximum entropy was found to be ",max_shared_binary_entropy," with frame size: ",optimal_frame_size
#     print "Now will be framming data with size: ", optimal_frame_size

    while alice_thread.event.is_set() or bob_thread.event.is_set():
        pass
    a_dict = alice_thread.binary_string_laser
    b_dict = bob_thread.binary_string_laser
    print "coinc window radius in bins", int(coincidence_window_radius/resolution)
    print "sync period size in bins ", int(sync_period/resolution)
    bob_full_dict = bob_corrected(bob[bob_chan == 4], sync_period, resolution, int(coincidence_window_radius/resolution))
#    
    print "before correction",len(intersect1d(bob[bob_chan == 4], alice[alice_chan == 0]))
    
    bob_corrected = do_correction(bob_full_dict, b_dict, a_dict, int(coincidence_window_radius/resolution))
#     print "bob_corrected"
    
#     savetxt("./DarpaQKD/bobCorrectedSystem.txt",json.dumps(bob_corrected))
    target = open("./DarpaQKD/bobCorrectedSystem.txt", 'w')
    target.write(json.dumps(a_dict))
    
    corrected_ttags = zeros(len(bob_corrected.keys()), dtype = uint64)
    i=0
    sync_block_size = int(sync_period/resolution)

    
    for key in bob_corrected.keys():
        corrected_ttags[i] = int(str(key))*100+int(float(str(bob_corrected[key])))
        i+=1
    indexes_of_order = corrected_ttags.argsort(kind = "mergesort")
#     self.channels = take(self.channels,indexes_of_order)
    corrected_ttags = take(corrected_ttags,indexes_of_order)
#     print corrected_ttags

    
    save("./DarpaQKD/bobNotCorrectedT.npy",bob)
    save("./DarpaQKD/bobCorrectedT.npy",corrected_ttags)
    save("./DarpaQKD/bobCorrectedC.npy",bob_chan)
    
    save("./DarpaQKD/aliceCorrectedT.npy",alice)
    save("./DarpaQKD/aliceCorrectedC.npy",alice_chan)

    
#     print "-->",corrected_ttags
    print "coincidences with correction !!!! ->",len(intersect1d(corrected_ttags, alice))
#     print get_coinc(alice, alice_chan, bob, bob_chan, 0, 4, int(coincidence_window_radius/resolution))

#     print "->>Coincidences", len(intersect1d(bob_corrected.values(), alice))/len(bob)
#     print bob_thread.binary_string_laser
#     print alice_thread.binary_string_laser
#     print "Making datasets equal..."
#     (alice_thread.frame_occupancies,bob_thread.frame_occupancies) = make_data_string_same_size(alice_thread.frame_occupancies,bob_thread.frame_occupancies)
#     (alice_thread.frame_locations,bob_thread.frame_locations) = make_data_string_same_size(alice_thread.frame_locations,bob_thread.frame_locations)
    alice_thread.do_set()
    bob_thread.do_set()
#     while alice_thread.event.is_set() or bob_thread.event.is_set():
#         pass
    
#     calculateStatistics(alice_thread.ttags,bob_thread.ttags,alice_thread.channels,bob_thread.channels, alice_thread.coincidence_window_radius,alice_thread.resolution)
    sys.stdout.flush()
#      
# #===================DEALS WITH OCCUPANCIES == 1==================================================

    print sum(intersect1d(alice_thread.binary_laser_string, bob_thread.binary_laser_string))
    mutual_frames_with_occupancy_one = logical_and(alice_thread.frame_occupancies==1,bob_thread.frame_occupancies==1)
#  
    alice_non_zero_positions_in_frame = alice_thread.frame_locations[mutual_frames_with_occupancy_one]
    bob_non_zero_positions_in_frame   = bob_thread.frame_locations[mutual_frames_with_occupancy_one]
#  
    print "Alice and Bob frame location DATA saved for LDPC procedure."
# #     fmt="%i" saves signed decimal integers
    savetxt("./resultsLaurynas/ALICE_BOB_NON_ZERO_POSITIONS_IN_FRAME1.csv",(alice_non_zero_positions_in_frame,bob_non_zero_positions_in_frame),fmt="%i")
    alice_thread.non_zero_positions = alice_non_zero_positions_in_frame
    bob_thread.non_zero_positions = bob_non_zero_positions_in_frame
    print "NONZERO: ", sum(bob_thread.non_zero_positions == alice_thread.non_zero_positions)," out of ", len(bob_thread.non_zero_positions)," % ", float(sum(bob_thread.non_zero_positions == alice_thread.non_zero_positions))/len(bob_thread.non_zero_positions)
#  =======================Will be announcing some part of the string==================================
    announce_fraction = 0.3
    print "Alice and Bob are now ANNOUNCING "+str(announce_fraction)+ " of their frame position strings"
    alice_thread.received_string = bob_thread.non_zero_positions[:int(len(bob_thread.non_zero_positions)*announce_fraction)]
    bob_thread.received_string = alice_thread.non_zero_positions[:int(len(alice_thread.non_zero_positions)*announce_fraction)]
    print "Succesfully ANNOUNCED will release threads", len(alice_thread.non_zero_positions),len(alice_thread.received_string)
    main_event.clear()
    
    LDPC_encode(alice_thread)
    
#=============Sending syndrome values and parity check matrix?=====
    bob_thread.syndromes = alice_thread.syndromes
    bob_thread.parity_matrix = alice_thread.parity_matrix
#==================================================================
    
    LDPC_decode(bob_thread)
        
#     stop = timeit.default_timer()
#     print stop - start 
      
# 
#     
#     
