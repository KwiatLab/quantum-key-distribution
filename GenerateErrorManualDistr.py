'''
Created on Aug 12, 2016

@author: laurynas
'''

from numpy import *
from random import randint,sample,shuffle

if __name__ == '__main__':
    dir = "./DarpaQKD/LDPC_alice_ttags8_1000.txt"
    alice = loadtxt(dir, dtype(uint8))
    bob = alice.copy()
    fraction_of_errors = 1.0
    alph = 8
    
    error_indices = zeros(len(alice),dtype = bool)
#     print int(fraction_of_errors*len(alice))
    error_indices[:int(fraction_of_errors*len(alice))] = True
    random.choice(arange(-1,2), p= [0.15,0.7,0.15])
    shuffle(error_indices)
    for i in range(len(error_indices)):
        if error_indices[i]:

            error = bob[i] + random.choice(arange(-1,2), p= [0.5,0.0,0.5])
            while (error >= alph or error < 0):
                error = bob[i] + random.choice(arange(-1,2), p= [0.5,0.0,0.5])
            bob[i] = error
#     print bob
#     for i in range(len(alice)):
#         print alice[i],"-",bob[i]
    savetxt("./DarpaQKD/LDPC_bob_ttags8_1000.txt",bob,fmt='%10d')   
    