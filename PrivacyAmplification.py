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
    
    for i in range (HASHED_KEY_LENGTH):
        xor_value = rand_matrix[i][0] ^ key[0] 
        for j in range (1, ORIGINAL_KEY_LENGTH):
            xor_value ^= rand_matrix[i][j] ^ key[j]
        hashed_key[i] = xor_value
    
    return hashed_key 