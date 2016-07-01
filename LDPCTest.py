'''
Created on Jul 1, 2016

@author: laurynas
'''
from SlepianWolf import encode,SW_LDPC
from SW_prep import transitionMatrix_data2,sequenceProbMatrix
from numpy import loadtxt
from ParityCheckMatrixGen import gallager_matrix


if __name__ == '__main__':
    alice = loadtxt("./DarpaQKD/LDPC_alice_ttags8.txt")
    bob = loadtxt("./DarpaQKD/LDPC_bob_ttags8.txt")
    
#  ========== encode ===================
    column_weight = 3
    number_parity_edges = 6
    frame_size = 8
    
    total_string_length = len(alice)
    
    number_of_parity_check_eqns_gallager = int(total_string_length*column_weight/number_parity_edges)
    parity_matrix = gallager_matrix(number_of_parity_check_eqns_gallager, total_string_length, column_weight, number_parity_edges)
    
    syndromes=encode(parity_matrix,alice,frame_size)
    print "syndromes: ", syndromes
#  ======================================
    
# ============ decode ===================
    decoder='bp-fft'
    iterations=70
    frozenFor=5
    
    
    transition_matrix = transitionMatrix_data2(bob,alice,frame_size)
    prior_probability_matrix = sequenceProbMatrix(alice,transition_matrix)
    belief_propagation_system = SW_LDPC(parity_matrix, syndromes, prior_probability_matrix, original=alice,decoder=decoder)
    
    belief_propagation_system.decode(iterations=iterations,frozenFor=frozenFor)
    
    
  