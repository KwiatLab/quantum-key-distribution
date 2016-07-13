'''
Created on Jul 12, 2016

@author: laurynas
'''
from numpy import *
if __name__ == '__main__':
    aliceCh_s = load("./Debugging/aliceCh.npy")
    aliceTtag_s = load("./Debugging/aliceTtags.npy")
    
    bobCh_s = load("./Debugging/bobCh.npy")
    bobTtag_s = load("./Debugging/bobTtags.npy")
    
    aliceCh_d = load("./Debugging/aliceChD.npy")
    aliceTtag_d = load("./Debugging/aliceTtagsD.npy")
    
    bobCh_d = load("./Debugging/bobChD.npy")
    bobTtag_d = load("./Debugging/bobTtagsD.npy")
    
#     print aliceCh_d == bobCh_d
    print (bobCh_d != bobCh_s)
    print aliceCh_d