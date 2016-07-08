'''
Created on Jul 8, 2016

@author: laurynas
'''
from numpy import *

if __name__ == '__main__':
    alice_ttags = load("./DarpaQKD/aliceTtagsBright.npy")
    alice_channels = load("./DarpaQKD/aliceChannelsBright.npy")
    
    bob_ttags = load("./DarpaQKD/bobTtagsBright.npy")
    bob_channels = load("./DarpaQKD/bobChannelsBright.npy")
    
    alice_ttags = alice_ttags[:len(alice_ttags)/10]
    alice_channels = alice_channels[:len(alice_channels)/10]

    
    bob_ttags = bob_ttags[:len(bob_ttags)/10]
    bob_channels = bob_channels[:len(bob_channels)/10]
    
    save("./DarpaQKD/aliceChannelsBright.npy",alice_channels)
    save("./DarpaQKD/aliceTtagsBright.npy",alice_ttags)
    
    save("./DarpaQKD/bobChannelsBright.npy",bob_channels)
    save("./DarpaQKD/bobTtagsBright.npy",bob_ttags)

    
    