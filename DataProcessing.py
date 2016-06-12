'''
Created on Feb 17, 2016

@author: laurynas
'''

import ttag
import numpy
from numpy import *
import sys


def read_raw_file(alice_raw_filename, name, resolution=None):
    alice_file_contents = loadtxt(alice_raw_filename)
    
    indexes_of_order = alice_file_contents[1, :].argsort(kind="mergesort")
    alice_file_contents = take(alice_file_contents, indexes_of_order, 1)
    
        
    alice__channels = alice_file_contents[0, :].astype(uint8);

    if (resolution):
#         print (alice_file_contents[1, :] - alice_file_contents[1, 0])/resolution
        alice_timetags = around((alice_file_contents[1, :] -
                                 alice_file_contents[1, 0]) / resolution).astype(uint64)
    else:
        alice_timetags = (alice_file_contents[1, :] -
                                alice_file_contents[1, 0]).astype(uint64)
#     print alice_timetags
    return (alice__channels,alice_timetags)
def read_processed_file(alice_processed_filename):

    alice_raw = load("./results/" + alice_processed_filename + ".npy")
    return (alice_raw)

def saveprep(fname,alice,alice_pol):
    print("Saving Data Arrays")
    sys.stdout.flush()
    save("./resultsLaurynas/"+name+"Ttags_"+fname+".npy",alice)
    save("./resultsLaurynas/"+name+"Channels_"+fname+".npy",alice_pol)
    
if (__name__ == '__main__'):
    
    #---------FOR LATER DEVELOPMENT-----------------------------
    # alice_raw_filename = raw_input("Insert filename with Alice counts: ")
    # bob_raw_filename = raw_input("Insert filename with Bob counts: ")
    numpy.set_printoptions(edgeitems = 100)
    #-----------------------------------------------------------
    
    
    
    alice_raw_filename = sys.argv[1]
    name = sys.argv[2]
    resolution = 156.25e-12
    
    alice_bob_ch_ttag =read_raw_file(alice_raw_filename,name ,resolution)
    saveprep("06032014_maxpower",alice_bob_ch_ttag[1],alice_bob_ch_ttag[0])
    
    print("Arrays are saved to file.")
    
    
