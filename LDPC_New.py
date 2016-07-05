'''
Created on Jul 5, 2016

@author: laurynas
'''
from operator import xor
from math import exp,log
from numpy import zeros,loadtxt, where,array,append
from ParityCheckMatrixGen import gallager_matrix
from SlepianWolf import encode

class CS_list_element(object):
    def __init__(self,bit_number, register_index):
        self.bit_number = bit_number
        self.LL_ptr = register_index


def LL_init():
    f_init = f_init_list[QBER]
    a_init = a_init_list[QBER]
    
    for i in range(1,k+1):
        LL_reg[i] = a_init

def a2f():
    table = zeros([Ma+1])
    for i in range(1,Ma+1):
        a = exp(-i/Na)
        f = (1+a)/(1-a)
        z = log(f)
        j = int(Nf*z+0.5)
        if j <= 0:
            j = 1
        if j > Mf:
            j = Mf
        table[i] = j

def f2a():
    table = zeros([Mf+1])
    for i in range(1,Mf+1):
        f = exp(i/Nf)
        a = (f-1)/(f+1)
        z = -log(a)
        j = int(Na*z+0.5)
        if j <= 0:
            j = 1
        if j > Ma:
            j = Ma
        f2a[i] = j
        

def cs_msgs_2_bits():
    m = #number of checksums
    
    for i in range(1,m+1):
        sign = d[i]
        big_alpha = 0
        j1 = cs_index[i]
        j2 = cs_index[i+1]
        
        for j in range(j1,j2):
            a1 = LL_reg[cs_list[j].LL_ptr]
            if (a1 < 0):
                sign = 1-sign
                big_alpha = big_alpha - a1
            else:
                big_alpha = big_alpha + a1
        
        for j in range(j1,j2):
            a1 = LL_reg[cs_list[j].LL_ptr]
            if (a1 < 0):
                p_sign = 1 - sign
                p_alpha = big_alpha + a1
            else:
                p_sign = sign
                p_alpha = big_alpha - a1
            
            if p_alpha <= 0:
                p_alpha = 1
            if p_alpha > Ma:
                p_alpha = Ma
            if p_sign == 0:
                LL_reg[cs_list[j].LL_ptr] = p_alpha
            else:
                LL_reg[cs_list[j].LL_ptr] = -p_alpha
                
def bit_msgs_2_cs():
    
    for i in range(n):
        j1 = LL_index[i]
        j2 = LL_index[i+1]
        f_tot = f_init
        
        for j in range(j1,j2):
            u = LL_reg[j]
            if (u > 0):
                u = a2f[u]
            else:
                u = -a2uf[-u]
            LL_reg[j] = u
            f_tot = f_tot + u
    
        for j in range(j1,j2):
            k = f_tot - LL_reg[j]
            if k < 0:
                p_sign =1
                k = -k
            else:
                p_sign = 0
            if k < 1:
                k = 1
            if k > Mf:
                k = Mf
            if p_sign == 1:
                LL_reg[j] = -f2a[k]
            else:
                LL_reg[j] = f2a[k]
        
        if f_tot < 0:
            y1[i] = 1 - y[i]
        else:
            y1[i] = y[i]
    
def converged():
    succes = 1
    for i in range(1,m+1):
        sum = 0
        j1 = cs_index[i]
        j2 = cs_index[i+1]
        for j in range(j1,j2):
            sum = xor(sum, y1[cs_list[j].bit_number])
        if sum != c[i]:
            success = 0
    
    return success

if __name__ == '__main__':
    LL_init()
    success = 0
    i = 0
    max_loops = 70
    #===========PROFILER 3 PARAMS ===========================
    
    Na = 555
    Ma = 3893
    Nf = 556
    Mf = 3896
    
    #========================================================
    alice = loadtxt("./DarpaQKD/LDPC_alice_ttags8.txt")
    y = loadtxt("./DarpaQKD/LDPC_bob_ttags8.txt")
    y1 = zeros(len(y))
    frame_size = 16
    number_of_parity_check_eqns_gallager = 6
    n = total_string_length = len(alice)
    column_weight = 3
    row_weight = 10
    number_parity_edges = 6
    
    parity_matrix = gallager_matrix(number_of_parity_check_eqns_gallager, total_string_length, column_weight, number_parity_edges)
    
    c = syndromes=encode(parity_matrix,alice,frame_size)
    b = encode(parity_matrix,y,frame_size)
    d = c^b
    m = len(c)
    
    
    
    #k = column_weight/m
    
    
    for j in range(1,m+1):
        cs_list = array(m,dtype = CS_list_element)
        cs_index = array(m)
        
        checksum_row = parity_matrix[j,:]
        bit_numbers = where(checksum_row == 1)
        
        for k in range(len(bit_numbers)):
            cs_list_el = CS_list_element(bit_numbers[k], register_value)
            cs_list[j] = append(cs_list[j], cs_list_el)
            
        cs_index[j] = cs_list[j]
    
    for i in range(max_loops):
        cs_msgs_2_bits()
        bit_msgs_2_cs()
        success = converged()
        i+=1        
        if success == 1:
            print "LDPC converged in",i,"loops\n"
            break
    
    if (success == 0):
        print "LDPC failed\n"
        