'''
Created on Jul 24, 2016

@author: laurynas
'''
from numpy import *
if __name__ == '__main__':
    at_data = load("./DarpaQKD/aliceTtagsBrightAttempt10th1.npy")
    bt_data = load("./DarpaQKD/bobTtagsBrightAttempt10th1.npy")
    
    ac_data = load ("./aliceChannelsBrightAttempt10th1.npy")
    bc_data = load("./bobChannelsBrightAttempt10th1.npy")
    
    laser_pulse = 260e-12
    resolution = 78.125e-12


    print("Saving Data Arrays")
    
    save("./DarpaQKD/aliceChannelsBrightAttempt10th1REBINNED.npy",ac_data)
    save("./DarpaQKD/aliceTtagsBrightAttempt10th1REBINNED.npy",(at_data*resolution/laser_pulse).astype(uint64))
    
    save("./DarpaQKD/bobChannelsBrightAttempt10th1REBINNED.npy",bc_data)
    save("./DarpaQKD/bobTtagsBrightAttempt10th1REBINNED.npy",(bt_data*resolution/laser_pulse).astype(uint64))