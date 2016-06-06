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
        self.event.set()
    def do_clear(self):
        self.event.clear()
    def do_set(self):
        self.event.set()
    def wait(self):
        self.event.wait()
    def run(self):
        while self.running:
            print self.name+" : Reading.csv files and converting to .npy\n"
            system("python ./DataProcessing.py "+self.raw_file_name+" "+self.name)
            
            print self.name+": Loading .npy data"
            (self.ttags,self.channels) = loadprep("06032014_maxpower",self.name)
            print "Loaded .npy data."
            
#           TODO: Add delays here
            
            print self.name +" will notify main of finish\n"
            self.do_clear()
            
            buf_num = ttag.getfreebuffer() 
            self.buffer = ttag.TTBuffer(buf_num,create=True,datapoints = len(self.ttags))
            self.buffer.resolution = self.resolution
            self.buffer.channels = max(self.channels)+1
 
         
            indexes_of_order = self.ttags.argsort(kind = "mergesort")
            self.channels = take(self.channels,indexes_of_order)
            self.ttags = take(self.ttags,indexes_of_order)
            
            #-----------------TODO: DEBUG ADD ARRAY-------------------
            # print("Alice ready. Adding Alice Data to Buffer")
#             print self.name+"Will try to add array to buffer"
#             self.buffer.addarray(self.channels,self.ttags)

            print "will be calculating optimal frame size but first lets notify that ready..." 
            while main_event.is_set():
                print "Waiting 1 "
            self.binary_string_laser = create_binary_string_from_laser_pulses(self.ttags)
            print self.name+" Calculating frame occupancies and locations:\n"
            self.frame_occupancies = calculate_frame_occupancy(self.binary_string_laser,self.frame_size)
            self.frame_locations = calculate_frame_locations_daniels_mapping(self.binary_string_laser, self.frame_occupancies, self.frame_size)
            self.event.clear()
            main_event.set()
            while main_event.is_set():
                print "Waiting 2"
            print self.name + " Finished."
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

    print "main is not blocked"
    while alice_thread.event.is_set() or bob_thread.event.is_set():
        main_event.wait()

    print "Making datasets equal"
    make_equal_size(alice_thread, bob_thread)

    alice_thread.do_set()
    bob_thread.do_set()
#     while alice_thread.event.is_set() or bob_thread.event.is_set():
#         main_event.wait()
#     alice_thread.wait()
#     bob_thread.wait()
#         
    print "STATISTICS: "
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
#     print "The maximum entropy was found to be ",max_shared_binary_entropy," with frame size: ",optimal_frame_size
#     
#     print "Now will be framming data with size: ", optimal_frame_size
    while alice_thread.event.is_set() or bob_thread.event.is_set():
        print "MAIN: wait till I can make datsets equal"
        main_event.wait()
    
    print "Making datasets equal..."
    (alice_thread.frame_occupancies,bob_thread.frame_occupancies) = make_data_string_same_size(alice_thread.frame_occupancies,bob_thread.frame_occupancies)
    (alice_thread.frame_locations,bob_thread.frame_locations) = make_data_string_same_size(alice_thread.frame_locations,bob_thread.frame_locations)
    print "Made datasets equal, RELEASING threads again!\n"
    main_event.clear()
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
    savetxt("./resultsLaurynas/ALICE_BOB_NON_ZERO_POSITIONS_IN_FRAME.csv",(alice_non_zero_positions_in_frame,bob_non_zero_positions_in_frame),fmt="%i")
    alice_thread.non_zero_positions = alice_non_zero_positions_in_frame
    bob_thread.non_zero_positions = bob_non_zero_positions_in_frame
 
#     stop = timeit.default_timer()
#     print stop - start 
      
# 
#     
#     