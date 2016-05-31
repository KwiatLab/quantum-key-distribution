import warnings
import numpy as np


def generate_matrix (rows, columns):
    return np.random.randint(2, size = (rows, columns))

def privacy_amplification (key, HASHED_KEY_LENGTH):
    
    ORIGINAL_KEY_LENGTH = len(key)
    hashed_key = np.zeros(HASHED_KEY_LENGTH, dtype=np.int)
    
    if (HASHED_KEY_LENGTH > len(key)):
        warnings.warn("the hashed key length is greater than original key")
    
    rand_matrix = generate_matrix(HASHED_KEY_LENGTH, ORIGINAL_KEY_LENGTH)
    
    for number_of_parity_check_eqns in range (HASHED_KEY_LENGTH):
        xor_value = rand_matrix[number_of_parity_check_eqns][0] ^ key[0] 
        for j in range (1, ORIGINAL_KEY_LENGTH):
            xor_value ^= rand_matrix[number_of_parity_check_eqns][j] ^ key[j]
        hashed_key[number_of_parity_check_eqns] = xor_value
    
    return hashed_key 