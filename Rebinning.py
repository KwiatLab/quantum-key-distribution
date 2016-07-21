'''
Created on Jul 20, 2016

@author: laurynas
'''
from numpy import *
from System import *

def rebinning(timetags,laser_pulse,resolution):
    return (timetags*resolution/laser_pulse).astype(uint64)
    
if __name__ == '__main__':
    channels1 = [0,1,2,3]
    channels2 = [4,5,6,7]
    laser_pulse = 260e-12
    resolution = 78.125e-12
    
    (alice,alice_ch) = load_data("alice",channels1,1000)
    (bob,bob_ch) = load_data("bob",channels2,1000)
    
    alice_new = rebinning(alice,laser_pulse,resolution)
    bob_nw = rebinning(bob, laser_pulse, resolution)
    
    print alice,"\n",alice_new
    
    