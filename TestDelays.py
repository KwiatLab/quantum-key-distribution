'''
Created on Jun 11, 2016

@author: laurynas
'''
from Delays import *
from numpy import *
if __name__ == '__main__':
    alice = loadtxt("./DataFiles/FakeDelays/AliceFake.txt")
    bob = loadtxt("./DataFiles/FakeDelays/BobFake.txt")
    
    aliceTtags= alice[1].astype(uint64)
    bobTtags=bob[1].astype(uint64)
    aliceChannels=alice[0].astype(uint8)
    bobChannels = bob[0].astype(uint8)
    
    
    resolution= 156.25e-12,
    channels1=channels2 = [2,3,4,5]

    (d1,d2) = calculate_delays(aliceTtags, aliceChannels, bobTtags, bobChannels,resolution= 156.25e-12,coincidence_window_radius = 5e-9, channel_array=channels1) 
    print d1,d2