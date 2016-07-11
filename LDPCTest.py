'''
Created on Jul 1, 2016

@author: laurynas
'''
from SlepianWolf import encode,SW_LDPC
from SW_prep import transitionMatrix_data,transitionMatrix_data2,transitionMatrix_data2_python,sequenceProbMatrix,randomMatrix,set_printoptions,nan
from numpy import loadtxt
from ParityCheckMatrixGen import gallager_matrix


if __name__ == '__main__':
    alice = loadtxt("./DarpaQKD/LDPC_alice_ttags8.txt")
    bob = loadtxt("./DarpaQKD/LDPC_bob_ttags8.txt")
    
#  ========== encode ===================
    column_weight = 2
    row_weight = 4
    frame_size = 16
    
    total_string_length = len(alice)
    print alice
    number_of_parity_check_eqns_gallager = int(total_string_length*column_weight/row_weight)
    parity_matrix = gallager_matrix(number_of_parity_check_eqns_gallager, total_string_length, column_weight, row_weight)
#     print parity_matrix.shape
#     parity_matrix = randomMatrix(total_string_length, number_of_parity_check_eqns_gallager, 4)
#     print parity_matrix
    syndromes=encode(parity_matrix,alice,alphabet=2)
#     print "syndromes: ", syndromes
#  ======================================
    
# ============ decode ===================
    decoder='log-bp-fft'
    iterations=70
    frozenFor=5
    
    set_printoptions(threshold=nan)
    transition_matrix = transitionMatrix_data2_python(alice,bob,2)
    prior_probability_matrix = sequenceProbMatrix(alice,transition_matrix)
#   
    print prior_probability_matrix.shape
    print "Transition matrix\n",transition_matrix
    belief_propagation_system = SW_LDPC(parity_matrix, syndromes, prior_probability_matrix, original=alice,decoder=decoder)
     
    print "KEY CORRECTNESS",sum(belief_propagation_system.decode(iterations=iterations,frozenFor=frozenFor) == alice)/float(len(alice))
    
    
  