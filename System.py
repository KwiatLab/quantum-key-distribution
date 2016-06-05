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
        self.event.set()
    def notify_finished(self):
        self.event.clear()
    def notify_ready(self):
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
            print self.name +" will notify main of finish\n"
            self.notify_finished()
            
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
            self.notify_finished()
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
# 
#     print "Alice: Reading.csv files and converting to .npy"
#     system("python ./DataProcessing.py "+alice_raw_filename+" alice")
#       
#     print "Alice: Loading .npy data"
#     (ttags,channels) = loadprep("06032014_maxpower","alice")
#      
#     
#     print "Bob: Reading.csv files and converting to .npy"
#     system("python ./DataProcessing.py "+bob_raw_filename+" bob")
#      
#     print "Bob: Loading .npy data"
#     (ttags,channels) = loadprep("06032014_maxpower","bob")
#     
   
    main_event = threading.Event()
    main_event.set()
    alice_thread.start()
    bob_thread.start()
#     event = threading.Event()
# #     event.set()
#     
#     print "joining threadds"
#     alice_thread.join()
#     bob_thread.join()
    print "main is not blocked"
    while alice_thread.event.is_set() or bob_thread.event.is_set():
        main_event.wait()
#         print alice_thread.event.is_set()

    print "Making datasets equal"
    make_equal_size(alice_thread, bob_thread)

    alice_thread.notify_ready()
    bob_thread.notify_ready()
#     while alice_thread.event.is_set() or bob_thread.event.is_set():
#         main_event.wait()
#     alice_thread.wait()
#     bob_thread.wait()
#         
    print "STATISTICS: "
    (alice,bob,alice_chan,bob_chan) = (alice_thread.ttags, bob_thread.ttags, alice_thread.channels, bob_thread.channels)
    statistics = calculateStatistics(alice,bob,alice_chan,bob_chan)
    print statistics
# 
#     stop = timeit.default_timer()
#     print stop - start 
     

    
    