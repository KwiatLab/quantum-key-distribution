'''
Created on Jun 28, 2016

@author: laurynas
'''
from System import get_timetag_corrections, do_correction
from numpy import loadtxt,zeros,sort, intersect1d,uint64,savetxt

if __name__ == '__main__':
    alice_ttags = loadtxt("./DataFiles/aliceTtags.txt")
    bob_ttags = loadtxt("./DataFiles/bobTtags.txt")
    resolution = 78.125e-12
    sync_period = 7.8125e-9
    coincidence_window_radius = 15
    
    
    a_dict = get_timetag_corrections(alice_ttags, resolution, sync_period, coincidence_window_radius)
    b_dict = get_timetag_corrections(bob_ttags, resolution, sync_period, coincidence_window_radius)
    print "alice ttags",alice_ttags
    print "bob ttags",bob_ttags
    print "alice dict",a_dict
    print "bob dict",b_dict
    bob_ttag_dict = {}
    sync_block_size = int(sync_period/resolution)
    print "sync_block_size",sync_block_size
    
    for i in range(len(bob_ttags)):
        ith = (bob_ttags[i] % sync_block_size)
        bob_ttag_dict[int(bob_ttags[i]/sync_block_size)] = ith
        
#     print bob_ttag_dict
    corrected = do_correction(bob_ttag_dict, b_dict, a_dict, coincidence_window_radius)
    new = zeros(len(corrected.keys()), dtype = uint64)
    i=0
#     resolution = 78.125e-12
#     sync_period = 8e-9

    for key in corrected.keys():
#         print sync_block_size % 10
#         print str(key), str(corrected[key]).zfill(2+sync_block_size % 10)
        new[i] = int(str(key))*100+int(float(str(corrected[key])))
        i+=1
#     for key in corrected.keys():
#         new[i] = key*int(sync_period/resolution) + corrected[key]
#         i+=1
#     print new
    print sort(new)
    print len(intersect1d(new,alice_ttags))
#     print a_dict
#     print b_dict