'''
Created on Jul 7, 2016

@author: laurynas
'''
from numpy import *

if __name__ == '__main__':
    file_names = ["buffer_0_dataset_1_time_10.csv",
                  "buffer_0_dataset_2_time_10.csv",
                  "buffer_0_dataset_3_time_10.csv",
                  "buffer_0_dataset_4_time_10.csv",
                  "buffer_0_dataset_5_time_10.csv",
                  "buffer_0_dataset_6_time_10.csv",
                  "buffer_1_dataset_1_time_10.csv",
                  "buffer_1_dataset_2_time_10.csv",
                  "buffer_1_dataset_3_time_10.csv",
                  "buffer_1_dataset_4_time_10.csv",
                  "buffer_1_dataset_5_time_10.csv",
                  "buffer_1_dataset_6_time_10.csv",
                  "buffer_2_dataset_1_time_10.csv",
                  "buffer_2_dataset_2_time_10.csv",
                  "buffer_2_dataset_3_time_10.csv",
                  "buffer_2_dataset_4_time_10.csv",
                  "buffer_2_dataset_5_time_10.csv",
                  "buffer_2_dataset_6_time_10.csv",
                  ]
    chan_main_data = array([])
    ttag_main_data = array([])
    offset = 0
    for file in file_names:
        if ("buffer_0") in file:
            offset = -1
        elif ("buffer_1") in file:
            offset = 2
        elif ("buffer_2") in file:
            offset = 2
        data = load(file)
        data[:,0] += offset
        
        chan_main_data = append(chan_main_data,data[:,0])
        ttag_main_data = append(ttag_main_data,data[:,1])
        
    indexes_of_order = ttag_main_data.argsort(kind = "quicksort")
    chan_main_data = take(chan_main_data,indexes_of_order)
    ttag_main_data = take(ttag_main_data,indexes_of_order)
    
    alice_channels = [0,1,2,3]
    bob_channels =   [4,5,6,7]
    
    save("./DarpaQKD/aliceChannelsBright.npy",chan_main_data[in1d(chan_main_data, alice_channels)])
    save("./DarpaQKD/aliceTtagsBright.npy",ttag_main_data[in1d(chan_main_data, alice_channels)])
    
    save("./DarpaQKD/bobChannelsBright.npy",chan_main_data[in1d(chan_main_data, bob_channels)])
    save("./DarpaQKD/bobTtagsBright.npy",ttag_main_data[in1d(chan_main_data, bob_channels)])
    
    
    