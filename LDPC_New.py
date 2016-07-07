'''
Created on Jul 5, 2016

@author: laurynas
'''
from operator import xor
from math import exp,log
from numpy import zeros,loadtxt, where,array,empty,append,int64, ndarray,uint64,sum
from ParityCheckMatrixGen import gallager_matrix
from SlepianWolf import encode

class CS_list_element(object):
    def __init__(self,bit_number, register_index):
        self.bit_number = bit_number
        self.LL_ptr = register_index


def LL_init():
    f_init = int(f_init_list[QBER], 0)
    a_init = int(a_init_list[QBER], 0)
    
    for i in range(k):
        LL_reg[i] = a_init
    
    return f_init,a_init,LL_reg

def get_a2f():
    table = zeros([Ma+1])
    for i in range(1,Ma+1):
        a = exp(-i/float(Na))
        f = (1+a)/(1-a)
        z = log(f)
        j = int(Nf*z+0.5)
        if j <= 0:
            j = 1
        if j > Mf:
            j = Mf
        table[i] = j        
    return table

def get_f2a():
    table = zeros([Mf+1])
    for i in range(1,Mf+1):
        f = exp(i/float(Nf))
#         print i/Nf
        a = (f-1)/(f+1)
#         print a
        z = -log(a)
        j = int(Na*z+0.5)
        if j <= 0:
            j = 1
        if j > Ma:
            j = Ma
        table[i] = j
    return table    

def cs_msgs_2_bits():
    
    for i in range(m):
        sign = d[i]
        big_alpha = 0
        
        j1 = cs_index[i]
            
        for j in range(len(j1)):
            a1 = LL_reg[cs_list[j].LL_ptr]
            if (a1 < 0):
                sign = 1-sign
                big_alpha = big_alpha - a1
            else:
                big_alpha = big_alpha + a1
        
        for j in range(len(j1)):
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
                
#     print "_",LL_reg
                
def bit_msgs_2_cs():
    
    for i in range(n):
#         print LL_index
        j1 = LL_index[i]
        
        f_tot = f_init
        
        for j in range(len(j1)):
            u = LL_reg[j]
            if (u > 0):
                u = a2f[u]
#                 print u 
            else:
                u = -a2f[-u]
            LL_reg[j] = u
#             print f_tot,u
            f_tot = f_tot + u
#             print f_tot
    
        for j in range(len(j1)):
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
    print "----"
    for i in range(m):
        sum = 0
        j1 = cs_index[i]
#         print 
        for j in range(len(j1)):
#             print y1[cs_list[j].bit_number]
#             print (alice ==  y1)
            pass
#             print y1[cs_list[j].bit_number]
        if sum != c[i]:
            success = 0
    
    return success

if __name__ == '__main__':
    
    a_init_list = [
                    "0x0004", # 1%
                    "0x0007", # 2%
                    "0x000B", # 3%
                    "0x000F", # 4%
                    "0x0012", # 5%
                    "0x0016", # 6%
                    "0x001A" # 7%
                    "0x001F", # 8%
                    "0x0023", # 9%
                    "0x0027", # 10%
                    "0x002B", # 11%
                    "0x0030" # 12%
                    ]
    
    f_init_list = [
                    "0x0316", # 1% QBER
                    "0x02A1", # 2%
                    "0x0259", # 3%
                    "0x0226", # 4%
                    "0x01FD", # 5%
                    "0x01DC", # 6%
                    "0x01C0", # 7%
                    "0x01A7", # 8%
                    "0x0190", # 9%
                    "0x017C", # 10%
                    "0x016A", # 11%
                    "0x0159" # 12%
                    ]
    #===========PROFILER 3 PARAMS ===========================
    
    Na = 555
    Ma = 3893
    Nf = 556
    Mf = 3896
    
    #========================================================
    
    
    alice = loadtxt("./DarpaQKD/LDPC_alice_ttags8.txt")
    y = loadtxt("./DarpaQKD/LDPC_bob_ttags8.txt")
    y1 = zeros(len(y), dtype = uint64)
    frame_size = 16
    n = total_string_length = len(alice)
    column_weight = 4
    row_weight = 8
    
    number_of_parity_check_eqns_gallager = int(total_string_length*column_weight/row_weight)
    
    parity_matrix = gallager_matrix(number_of_parity_check_eqns_gallager, total_string_length, column_weight, row_weight)
#     print parity_matrix
#     print parity_matrix
    c = syndromes=encode(parity_matrix,alice,frame_size)
    b = encode(parity_matrix,y,frame_size)
    d = c^b
    m = len(c)
    a2f = get_a2f()
    f2a = get_f2a()
    
    QBER = 10
    k = int(row_weight*m) #Bits per checksum
#     print k
#     print sum(sum(parity_matrix))
    LL_reg = empty(k, dtype = int64)

    (f_init, a_init, LL_reg) = LL_init()
    
    success = 0
    i = 0
    max_loops = 70

        
    cs_list = array([], dtype = CS_list_element)
    cs_index = []
    
    reg_index = 0
    for j in range(m):
        
        checksum_row = parity_matrix[j,:]
        bit_numbers = where(checksum_row == 1)[0]
        this_row_weight = len(bit_numbers)

        cs_index_el = array(this_row_weight, dtype = CS_list_element)
        for k in range(this_row_weight):
            cs_list_array_el = CS_list_element(bit_numbers[k], LL_reg[reg_index])
            reg_index +=1
            cs_index_el =append(cs_index_el,cs_list_array_el)

            cs_list = append(cs_list,cs_list_array_el)
        cs_index.append(cs_index_el)
#         print cs_index
 
    LL_index = []
    relative_index = 0
    for i in range(n):
        bit_column = parity_matrix[:,i]
        this_column_weight = len(where(bit_column == 0)[0])
        LL_index.append(LL_reg[relative_index:relative_index+this_column_weight]) 
        relative_index += this_column_weight

    for i in range(max_loops):
        cs_msgs_2_bits()
        bit_msgs_2_cs()
        print LL_reg
        success = converged()
        i+=1        
        if success == 1:
            print "LDPC converged in",i,"loops\n"
            break
    
    if (success == 0):
        print "LDPC failed\n"
        