'''
Created on May 17, 2016

@author: laurynas
'''
from entropy_calculator import *
from numpy import *

if __name__ == '__main__':
    binary_string_alice = array([4,7,8,12,15,18])
    binary_string_bob = array([4,8,13,15,19])
    frame_size = 16

    entropy_calc(binary_string_alice,binary_string_bob, frame_size)

    