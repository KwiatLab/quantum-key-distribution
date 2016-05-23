'''
Created on Feb 17, 2016

@author: laurynas
'''

import ttag
import numpy
from numpy import *
import sys


def read_raw_file(alice_raw_filename, bob_raw_filename, resolution=None):
    alice_file_contents = loadtxt(alice_raw_filename)
    bob_file_contents = loadtxt(bob_raw_filename)
    
    indexes_of_order = alice_file_contents[1, :].argsort(kind="mergesort")
    alice_file_contents = take(alice_file_contents, indexes_of_order, 1)
    indexes_of_order = bob_file_contents[1, :].argsort(kind="mergesort")
    bob_file_contents = take(bob_file_contents, indexes_of_order, 1)
    
        
    alice__channels = alice_file_contents[0, :].astype(uint8);
    bob__channels = bob_file_contents[0, :].astype(uint8);

    if (resolution):
        alice_timetags = around((alice_file_contents[1, :] -
                                 alice_file_contents[1, 0]) / resolution).astype(uint64)
        bob_timetags   = around((bob_file_contents[1, :] -
                                 bob_file_contents[1, 0]) / resolution).astype(uint64)
    else:
        alice_timetags = (alice_file_contents[1, :] -
                                alice_file_contents[1, 0]).astype(uint64)
        bob_timetags =   (bob_file_contents[1, :] -
                                bob_file_contents[1, 0]).astype(uint64)
    print("alice channels")
    print(alice__channels)
    print("bob channels")
    print(bob__channels)
    print("alice ttags")
    print(alice_timetags)
    print("bob ttags")
    print(bob_timetags)
    
    return ([alice__channels,alice_timetags,bob__channels,bob_timetags])
def read_processed_file(alice_processed_filename, bob_processed_filename):

    alice_raw = load("./results/" + alice_processed_filename + ".npy")
    bob_raw = load("./results/" + bob_processed_filename + ".npy")
    return (alice_raw, bob_raw)

def saveprep(fname,alice,bob,alice_pol,bob_pol):
    print("Saving Data Arrays")
    sys.stdout.flush()
    save("./resultsLaurynas/aliceTtags_"+fname+"_full.npy",alice)
    save("./resultsLaurynas/bobTtags_"+fname+"_full.npy",bob)
    save("./resultsLaurynas/aliceChannels_"+fname+"_full.npy",alice_pol)
    save("./resultsLaurynas/bobChannels_"+fname+"_full.npy",bob_pol)
    
if (__name__ == '__main__'):
    
    #---------FOR LATER DEVELOPMENT-----------------------------
    # alice_raw_filename = raw_input("Insert filename with Alice counts: ")
    # bob_raw_filename = raw_input("Insert filename with Bob counts: ")
    numpy.set_printoptions(edgeitems = 100)
    #-----------------------------------------------------------
    
    
    
    alice_raw_filename = "DataFiles/06032014_maxpower_268_0.csv"
    bob_raw_filename = "DataFiles/06032014_maxpower_268_1.csv"
    
    resolution = 5e-11 
    alice_bob_ch_ttag =read_raw_file(alice_raw_filename, bob_raw_filename,resolution)
    saveprep("06032014_maxpower",alice_bob_ch_ttag[1],alice_bob_ch_ttag[3],alice_bob_ch_ttag[0],alice_bob_ch_ttag[2])
    
    print("Arrays are saved to file.")
    
    
