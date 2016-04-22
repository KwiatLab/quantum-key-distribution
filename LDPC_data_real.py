from numpy import *
import numpy
from matplotlib import *
from multiprocessing import Pool
from scipy.optimize import curve_fit
import ttag
import sys
from scipy.weave import inline
#import graphs
from SW_prep import *
from SlepianWolf import *
from entropy_calculator import *



def createLDPCdata(timetags,polarizations,total_number_of_frames=None,frame_size=16):
    frame_numbers = timetags.copy()
    frame_numbers /= frame_size

    #Get the number of time bins
    if (total_number_of_frames == None): total_number_of_frames = frame_numbers[-1]+1

    #Create an array for the frame occupancy and location sequences
    frame_occupancy = zeros(total_number_of_frames,dtype=uint32)
    frame_location = zeros(total_number_of_frames,dtype=uint64)
    person_p = zeros(total_number_of_frames,dtype=uint8)
    frame_size = int(frame_size)

    print(frame_numbers)
    number_of_timetags = int(argmax(frame_numbers>total_number_of_frames))
    if (number_of_timetags==0 and frame_numbers[0] > total_number_of_frames):
        print "ERROR: TOT is smaller than minvalue!"
        return None
    if (number_of_timetags==0):
        number_of_timetags = int(len(timetags))
    
    print("number of timetags: ",number_of_timetags)
    maxtag = (timetags[number_of_timetags-1])
    print("total number of frames: ",total_number_of_frames)

    code = """
        long long z = 0;
        for (;z<number_of_timetags;z++) {
            frame_occupancy[frame_numbers[z]] +=1;
            frame_location[frame_numbers[z]] += (timetags[z] % frame_size)*(pow((double)frame_size,(double)(frame_occupancy[frame_numbers[z]]-1)));
            person_p[frame_numbers[z]] = polarizations[z];
        }
    """
    inline(code,["frame_occupancy","frame_location","person_p","frame_numbers","polarizations","timetags","number_of_timetags","frame_size"],headers=["<math.h>"])

    """
    for z in xrange(len(timetags)):
        frame_occupancy[frame_numbers[z]] += 1
        frame_location[frame_numbers[z]] = (timetags[z] % frame_size) * frame_occupancy[frame_numbers[z]]*frame_size
        person_p[frame_numbers[z]] = polarizations[z]
    """
    return (frame_occupancy,frame_location,person_p,maxtag)

