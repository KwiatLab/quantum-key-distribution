'''
Created on May 15, 2016

@author: laurynas
'''

import numpy
from numpy import *
import sys

def read_raw_file(alice_raw_filename, bob_raw_filename):
    print "setting options"
    numpy.set_printoptions(edgeitems = 100)
    print "reading alice..."
    alice_file_contents = loadtxt(alice_raw_filename+".csv")
    print "reading bob..."
    bob_file_contents = loadtxt(bob_raw_filename+".csv")
    print "printing alice..."
    print (alice_file_contents)
    print "printing done"
    return (alice_file_contents,bob_file_contents)

    # indexes_of_order = alice_file_contents[1, :].argsort(kind="mergesort")
    # alice_file_contents = take(alice_file_contents, indexes_of_order, 1)
    # indexes_of_order = bob_file_contents[1, :].argsort(kind="mergesort")
    # bob_file_contents = take(bob_file_contents, indexes_of_order, 1)
    
        
    # alice__channels = alice_file_contents[0, :].astype(uint8);
    # bob__channels = bob_file_contents[0, :].astype(uint8);

def trim_files(file_name, file_contents, trim_factor = 50):
    flag = True
    for element in file_contents:
#         print element
        new_size = len(element)/trim_factor
#         print new_size
        element_store = element[:new_size]
#         print element_store
        if flag:
            new_file_content = array(element_store)
            flag = False
        else:
            new_file_content = vstack((new_file_content,array(element_store)))
    savetxt(file_name + "_trimmed.csv", new_file_content)

if __name__ == '__main__':
    numpy.set_printoptions(edgeitems = 100)
    #-----------------------------------------------------------
    
    
     
    alice_raw_filename = "06032014_maxpower_268_0"
    bob_raw_filename = "06032014_maxpower_268_1"
    # print "will start to read"
    (alice_file_contents,bob_file_contents) = read_raw_file(alice_raw_filename, bob_raw_filename)
    # trims by default by factor of 50
    trim_files(alice_raw_filename,alice_file_contents)
    trim_files(bob_raw_filename, bob_file_contents)
    print "Printing trimmed alice:"
    print (loadtxt(alice_raw_filename+"_trimmed.csv"))
    print "Printing trimmed bob:"
    print (loadtxt(bob_raw_filename+"_trimmed.csv"))
    
#     print (loadtxt("alice_trimmed.csv"))