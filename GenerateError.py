'''
Created on Aug 12, 2016

@author: laurynas
'''

from numpy import *
from random import randint,sample,shuffle

if __name__ == '__main__':
    dir = "./DarpaQKD/LDPC_alice_ttags8_100.txt"
    alice = loadtxt(dir, dtype(uint8))
    bob = alice.copy()
    fraction_of_errors = 0.99
    alph = 8
    
    error_indices = zeros(len(alice),dtype = bool)
#     print int(fraction_of_errors*len(alice))
    error_indices[:int(fraction_of_errors*len(alice))] = True

    shuffle(error_indices)
    for i in range(len(error_indices)):
        if error_indices[i]:
            pos_poisson_string  = random.poisson(1,(1,1,alph))
            full_poisson_string = append(pos_poisson_string,pos_poisson_string*(-1))
            full_poisson_string[where(full_poisson_string == 0)] = sample(set([-1,1]),1)
            error = bob[i] + full_poisson_string[randint(0,len(full_poisson_string)-1)]
            while (error >= alph or error < 0):
                error = bob[i] + full_poisson_string[randint(0,len(full_poisson_string)-1)]
            bob[i] = error
            
    savetxt("./DarpaQKD/LDPC_bob_ttags8_100.txt",bob,fmt='%10d')   
    