'''
Created on Jul 11, 2016

@author: laurynas
'''
import pyldpc
from numpy import *

if __name__ == '__main__':
    msg = array([0,1,0,1,1,0])
    
    n = 15  # Number of columns
    column_weight = 4 # Number of ones per column, must be lower than d_c (because H must have more rows than columns)
    row_weight = 5 
    snr = 8
    
    H = pyldpc.RegularH(n,column_weight,row_weight)
    tG = pyldpc.CodingMatrix(H)
    
    new_H,sys_tG = pyldpc.CodingMatrix_systematic(H)
    
    y = pyldpc.Coding(tG,msg,snr)
    x_decoded = pyldpc.Decoding_logBP(H,y,snr,5)
    print("H.x' = ",pyldpc.BinaryProduct(H,x_decoded))
    
    v_received = pyldpc.DecodedMessage(tG,x_decoded)
    print v_received
    