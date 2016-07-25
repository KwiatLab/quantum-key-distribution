'''
Created on Jul 24, 2016

@author: laurynas
'''
from numpy import *
if __name__ == '__main__':
    data = loadtxt("./DarpaQKD/Alice1_Bob1.csv")
    alice_channels = [0,1,2,3]
    bob_channels = [4,5,6,7]
    laser_pulse = 260e-12
    resolution = 78.125e-12
    channels = data[:,0]
    timetags = data[:,1]
    timetags = (timetags*resolution/laser_pulse).astype(uint64)

    print("Saving Data Arrays")
    
    save("./DarpaQKD/aliceChannelsREBINNEDfull.npy",channels[in1d(channels, alice_channels)])
    save("./DarpaQKD/aliceTtagsREBINNEDfull.npy",timetags[in1d(channels, alice_channels)])
    
    save("./DarpaQKD/bobChannelsREBINNEDfull.npy",channels[in1d(channels, bob_channels)])
    save("./DarpaQKD/bobTtagsREBINNEDfull.npy",timetags[in1d(channels, bob_channels)])