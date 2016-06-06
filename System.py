'''
Created on Jun 4, 2016

@author: laurynas
'''
import sys
import threading
from os import system
import DataProcessing
from numpy import *
import ttag
from multiprocessing import *
from subprocess import Popen, list2cmdline
from Statistics import *

import timeit

def make_equal_size(alice_thread,bob_thread):
    if (len(alice_thread.ttags) > len(bob_thread.ttags)):
        alice_thread.ttags = alice_thread.ttags[:bob_thread.ttags]
        alice_thread.channels = alice_thread.channels[:len(bob_thread.channels)]
    else:
        bob_thread.ttags    = bob_thread.ttags[:len(alice_thread.ttags)]
        bob_thread.channels = bob_thread.channels[:len(alice_thread.channels)]
def loadprep(fname,name):

    sys.stdout.flush()
    alice = load("./resultsLaurynas/resultsLaurynas/"+name+"Ttags_"+fname+".npy")
    alice_pol=load("./resultsLaurynas/resultsLaurynas/"+name+"Channels_"+fname+".npy")
    return (alice,alice_pol)

class PartyThread(threading.Thread):
    def __init__(self, resolution, raw_filename,name):
        threading.Thread.__init__(self)
        self.running = True
        self.name = name
        self.resolution = resolution
        self.raw_file_name = raw_filename
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
        self.event.set()
    def do_clear(self):
        self.event.clear()
    def do_set(self):
        self.event.set()
    def run(self):
        while self.running:
            print self.name.upper()+" : Reading.csv files and converting to .npy\n"
            system("python ./DataProcessing.py "+self.raw_file_name+" "+self.name)
            print self.name.upper()+": Loading .npy data"
            (self.ttags,self.channels) = loadprep("06032014_maxpower",self.name)
#           TODO: Add delays here
            print self.name.upper() +" FINISHED with data. Will notify main.\n"        
            buf_num = ttag.getfreebuffer() 
            self.buffer = ttag.TTBuffer(buf_num,create=True,datapoints = len(self.ttags))
            self.buffer.resolution = self.resolution
            self.buffer.channels = max(self.channels)+1
 
         
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
            self.binary_string_laser = create_binary_string_from_laser_pulses(self.ttags)
            print self.name+" Calculating frame occupancies and locations:\n"
            self.frame_occupancies = calculate_frame_occupancy(self.binary_string_laser,self.frame_size)
            self.frame_locations = calculate_frame_locations_daniels_mapping(self.binary_string_laser, self.frame_occupancies, self.frame_size)
            print self.name.upper()+": Notifying main that finished with laser strings. Will be waiting for"

#===========READY TO ANNOUNCE======================
            self.do_clear()
            while main_event.is_set():
                pass
            print self.name.upper()+": I was released and will do error calc.\n"    
            self.error_rate = sum(self.received_string == self.non_zero_positions[:len(self.received_string)])/len(self.received_string)
            print self.name.upper() + " Finished with error rate: ",self.error_rate,"\n"
            
            
            self.running = False
               
if __name__ == '__main__':

    
    alice_raw_filename = "./DataFiles/ShorterFiles/06032014_maxpower_268_0_trimmed.csv"
    bob_raw_filename = "./DataFiles/ShorterFiles/06032014_maxpower_268_1_trimmed.csv"
    
    set_printoptions(edgeitems = 100)
    resolution = 5e-11
    
    alice_event = threading.Event()
    alice_event.set()
   
    bob_event = threading.Event()
    bob_event.set() 
    
    alice_thread = PartyThread(resolution, alice_raw_filename, "alice")
    bob_thread = PartyThread(resolution, bob_raw_filename,"bob")
    start = timeit.default_timer()  
   
    main_event = threading.Event()
    main_event.set()
    alice_thread.start()
    bob_thread.start()

    print "MAIN: will wait till AB finished loading data."
    while alice_thread.event.is_set() or bob_thread.event.is_set():
        pass
    
    print "MAIN: Making datasets equal."
    make_equal_size(alice_thread, bob_thread)
    
    print "MAIN: STATISTICS: "
    (alice,bob,alice_chan,bob_chan) = (alice_thread.ttags, bob_thread.ttags, alice_thread.channels, bob_thread.channels)
#     statistics = calculateStatistics(alice,bob,alice_chan,bob_chan)
#     print statistics
    
#     max_shared_binary_entropy = max(statistics.values())
#     optimal_frame_size = int(list(statistics.keys())[list(statistics.values()).index(max_shared_binary_entropy)])
    optimal_frame_size = 512
    alice_thread.frame_size = optimal_frame_size
    bob_thread.frame_size = optimal_frame_size
    
    print "MAIN: Optimal size calculated and set for both threads, release THEM!"
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
    
    
    print "Making datasets equal..."
    (alice_thread.frame_occupancies,bob_thread.frame_occupancies) = make_data_string_same_size(alice_thread.frame_occupancies,bob_thread.frame_occupancies)
    (alice_thread.frame_locations,bob_thread.frame_locations) = make_data_string_same_size(alice_thread.frame_locations,bob_thread.frame_locations)
    alice_thread.do_set()
    bob_thread.do_set()
#     while alice_thread.event.is_set() or bob_thread.event.is_set():
#         pass
    
    
    sys.stdout.flush()
#      
# #===================DEALS WITH OCCUPANCIES == 1==================================================
#  
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
#  =======================Will be announcing some part of the string==================================
    announce_fraction = 0.1
    print "Alice and Bob are now ANNOUNCING "+str(announce_fraction)+ " of their frame position strings"
    announce_fraction = 0.1
    alice_thread.received_string = bob_thread.non_zero_positions[:int(len(bob_thread.non_zero_positions)*announce_fraction)]
    bob_thread.received_string = alice_thread.non_zero_positions[:int(len(alice_thread.non_zero_positions)*announce_fraction)]
    print "Succesfully ANNOUNCED will release threads", len(alice_thread.non_zero_positions),len(alice_thread.received_string)
    main_event.clear()
        
#     stop = timeit.default_timer()
#     print stop - start 
      
# 
#     
#     
