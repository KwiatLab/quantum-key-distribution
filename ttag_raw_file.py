'''
Created on Jun 16, 2016

@author: laurynas
'''
/* The UIUC/NCSA license:

Copyright (c) 2014 Kwiat Quantum Information Group
All rights reserved.

Developed by:    Kwiat Quantum Information Group
University of Illinois, Urbana-Champaign (UIUC)
http://research.physics.illinois.edu/QI/Photonics/

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal with the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers
in the documentation and/or other materials provided with the distribution.

Neither the names of Kwiat Quantum Information Group, UIUC, nor the names of its contributors may be used to endorse
or promote products derived from this Software without specific prior written permission.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
*/

/**************************************************************************************

                            Time Tagger Data Analysis
                                  Implementation

**************************************************************************************/
/*
Note: Documentation in ttag.h
*/

#define TTAG_ENABLE_INTERNAL_MACROS //Enables macros used for simple output of warnings/errors/assertions
#include "ttag.h"
#include <stdlib.h>
#include <stdio.h>


/*
Internal functions: No checks for validity happen here, they are used for simplifying some of the analysis functions
*/
static inline uint64_t ttinternal_getNextChannelIndex(const tt_buf *const buffer, uint64_t curindex, const uint8_t channel, const uint64_t mintime, const uint64_t minindex) {
    if (curindex > minindex) {
        //Note the curindex > minindex -> it does not get to the last, since minindex can be 0
        for (--curindex; curindex > minindex && tt_tag(buffer,curindex) >= mintime; curindex--) {
            if (tt_channel(buffer,curindex)==channel) {
                //The minindex that was passed in can be invalid if data was acquired during this function's execution time.
                //    handle this eventuality here.
                if (curindex > tt_minindex(buffer)) {
                    return curindex;
                } else {
                    return ~((uint64_t)0);
                }
            }
        }
        //Now check for the last 0
        if (curindex == minindex && curindex >= tt_minindex(buffer) && tt_tag(buffer,curindex) >= mintime) {
            if (tt_channel(buffer,curindex)==channel) {
                return curindex;
            }
        }
    }
    return ~((uint64_t)0);
}

static inline int ttinternal_getMinChannel(const tt_buf *const buffer, const uint64_t *const channelMax, const uint64_t *const channelIndex, const uint8_t channels) {
    int i,mini;
    uint64_t minD = ~((uint64_t)0);
    mini = -1;
    for (i=0;i < channels;i++) {
        if (channelIndex[i]!= ~((uint64_t)0) && channelMax[i] - tt_tag(buffer,channelIndex[i]) < minD) {
            mini = i;
            minD= channelMax[i] - tt_tag(buffer,channelIndex[i]);
        }
    }
    return mini;
}

/*
TT_SINGLES
    Given a time, it returns the number of counts in the most recent time period. An optional singles
    array of size number-of-channels can be passed in to get the singles per channel in the time.
    The channel array is NOT zeroed inside the function, so it is your job to either zero it before
    passing it in, or doing whatever else you want. whatever contents of the singles array will be
    incremented for each count (a feature, not a bug :-P).
    
    The raw function gives singles for timebins back from datapoint at dataindex, while the tt_singles
    function gives the singles back from the reference time.
*/
TT_DEF_ uint64_t tt_rawsingles(const tt_buf *const buffer, uint64_t timebins, uint64_t* singlesarray, uint64_t dataindex) {
    uint64_t singles = 0;
    
    //TT_ASSERT(buffer,0);
    TT_ASSERT(timebins,0);
    TT_ASSERT(dataindex < tt_datanum(buffer),0);
    TT_ASSERT(dataindex >= tt_minindex(buffer),0);
    //If we are to populate a singles array for each channel, we do it now, knowing the total
    //  number of datapoints
    if (singlesarray) {
        //Here we check the singles and write a singles array:
        uint64_t i;
        
        //If the requested time is greater or equal to the total time, return all the buffered datapoints.
        TT_WARN(timebins+tt_tag(buffer,tt_minindex(buffer)) >= tt_tag(buffer,dataindex),,
            "Requested time is more than entire buffer, singles won't be correct!");
        
        for (i=dataindex; tt_tag(buffer,i)+timebins >= tt_tag(buffer,dataindex) && tt_tag(buffer,i) <= tt_tag(buffer,dataindex); i--) {
            singles++;
            singlesarray[tt_channel(buffer,i)]++;
            if (i<=tt_minindex(buffer)) break;
        }
        
    } else {
        //If the per-channel array is not needed, use a faster way to calculate the singles
        //The number of singles returned by bins2points is excluding the current datapoint. We want to include it.
        singles = tt_bins2points(buffer,dataindex,timebins)+1;
    }
    
    return singles;
}

TT_DEF_ uint64_t tt_singles(const tt_buf *const buffer, double time,uint64_t* singles) {
    //TT_ASSERT(buffer,0);
    TT_ASSERT(time,0);
    
    //Check for reference time
    TT_ASSERT(!isnan(tt_resolution(buffer)),0);
    
    //Make sure that there is data - and if not, return empty result
    TT_WARN(tt_datanum(buffer)==0,return 0;,"Buffer is empty!");

    return tt_rawsingles(buffer,tt_subtractreference(buffer,tt_time2bin(buffer,time)),singles,tt_datanum(buffer)-1);
}

/*
TT_COINCIDENCES
    Given a time, a radius within which to look, as well as an optional array of delays per channel, this function
    returns a matrix of coincidences of each channel with each channel, as well as an optional array of singles per
    channel within the period, and taking into account the delays between channels. tt_rawcoincidences also returns an
    error code
*/

TT_DEF_ uint64_t* tt_rawcoincidences(const tt_buf *const buffer, uint64_t timebins, uint64_t radius, uint64_t* coincidenceMatrix, uint64_t* delayArray, uint64_t dataindex) {
    uint64_t* channelIndex; 
    uint64_t* channelMin;
    uint64_t* channelMax;
    uint64_t icounter;
    int i,j;
    uint8_t channels;
    
    //Make sure everything is working correctly
    //TT_ASSERT(buffer,NULL);
    TT_ASSERT(timebins,NULL);
    TT_ASSERT(dataindex < tt_datanum(buffer),NULL);
    
    //Create simplifying shortcut to channels
    channels = tt_channels(buffer);
    
    
    //This all can be turned into one memory allocation, at the cost of readability
    channelIndex = (uint64_t*)malloc(sizeof(uint64_t)*channels);
    TT_CHKERR(channelIndex,return NULL,"Memory allocation failed!");
    channelMin = (uint64_t*)malloc(sizeof(uint64_t)*channels);
    TT_CHKERR(channelMin,free(channelIndex);return NULL,"Memory allocation failed!");
    channelMax = (uint64_t*)malloc(sizeof(uint64_t)*channels);
    TT_CHKERR(channelMax,free(channelMin); free(channelIndex); return NULL,"Memory allocation failed!");
    
    //Allocate a coincidence matrix if one is not given
    if (!coincidenceMatrix) {
        coincidenceMatrix = (uint64_t*)calloc(channels*channels,sizeof(uint64_t));
        
        TT_CHKERR(coincidenceMatrix, free(channelMax); free(channelMin); free(channelIndex); return NULL,"Failed to allocate coincidence matrix.");
    }
    
    
    
    //Find channel iterator locations, and set the bounds for channel times (channelMin, channelMax)
    for (i=0; i< channels;i++) {
    
        //put the index at current location
        channelIndex[i] = dataindex;
        
        //Set the channel's maximum timing
        channelMax[i] = tt_tag(buffer,dataindex);
        
        //If there is a delay array, find iterator locations offset by the delay amounts
        if (delayArray) {
            channelIndex[i] -= tt_bins2points(buffer,dataindex,delayArray[i]);
            //The difference between bins2points and the intended result is that they behave differently when there are multiple points with the same time tag
            //  on the boundary
            //For example, with boundary being 1200,  the datapoint returned is:
            //         100,1000,1200,1200,1200,1200,1400,1800
            //resut:            ^
            //want:                            ^
            //These are, of course, details, but they are important nevertheless! We have to account for them!
            for (icounter=channelIndex[i];icounter < dataindex && tt_tag(buffer,icounter+1)==tt_tag(buffer,channelIndex[i]);icounter++);
            channelIndex[i]=icounter;
            
            //Make sure that the channel's Max time does not overflow. Look up TT_WARN to find
            //  out why putting an else there is kosher.
            TT_WARN(channelMax[i] <= delayArray[i],channelMax[i]=0,"Delay for channel %i goes beyond available time!",i) 
                else {
                    channelMax[i] -= delayArray[i];
                }
        }
        
        //Next, subtract the allowed time to get the minimum time tag to accept. The minindex is taken into account later
        TT_WARN(channelMax[i] <= timebins, channelMin[i] = 0,"Time for channel %i goes beyond total data!",i)
            else {
                channelMin[i] = channelMax[i] - timebins;
            }
    }

    //Now make sure that the channelIndices are pointing to the correct channels, and that they are within their delays
    for (i=0;i< channels;i++) {
        
        //Make sure that the channel index is within channelMax
        while (((~((uint64_t)0))!=channelIndex[i]) && tt_tag(buffer,channelIndex[i]) > channelMax[i]) {
            channelIndex[i] = ttinternal_getNextChannelIndex(buffer,channelIndex[i], (uint8_t)i,channelMin[i],tt_minindex(buffer));
        }
        
        if (((~((uint64_t)0))!=channelIndex[i]) && i!=(int)tt_channel(buffer,channelIndex[i])) {
            //Not the correct channel, move index pointer to the correct one.
            
            channelIndex[i] = ttinternal_getNextChannelIndex(buffer,channelIndex[i], (uint8_t)i,channelMin[i],tt_minindex(buffer));
        }
    }
    //Preliminaries are set. This means:
    // - channelIndex contains the index of the first element of each channel
    // - channelMax contains the time from which to count down 'timebins' for each time
    // - channelMin is the minimum time tag to allow for each channel
    
    //Now, we take 
    while ((i= ttinternal_getMinChannel(buffer,channelMax,channelIndex,channels))!=-1) {
        
        //Check for coincidences. At each point in time, i is the index of the channel with the highest "time". We take i, and compare to each of the other channels,
        //  to see if they are within radius. For each that is within radius, we add a coincidence at both locations in the coincidence matrix.
        for (j=0; j< channels; j++) {            
            if (((~((uint64_t)0))!=channelIndex[j]) && (channelMax[i]-tt_tag(buffer,channelIndex[i])+radius >= channelMax[j]-tt_tag(buffer,channelIndex[j]))) {
                coincidenceMatrix[channels*i+j]++;
                if (i!=j) coincidenceMatrix[channels*j+i]++;    //Make sure we don't double-add on diagonal
            }
        }
        
        channelIndex[i] = ttinternal_getNextChannelIndex(buffer,channelIndex[i], (uint8_t)i, channelMin[i],tt_minindex(buffer));
    }
    
    //Free all necessary memory
    free(channelIndex);
    free(channelMin);
    free(channelMax);
    
    return coincidenceMatrix;
}

TT_DEF_ uint64_t* tt_coincidences(const tt_buf *const buffer, double time, double radius, uint64_t* coincidenceMatrix, double* delayArray) {
    uint64_t* delays = NULL;
    
    //TT_ASSERT(buffer,NULL);
    TT_ASSERT(time,NULL);
    
    //Check for reference time
    TT_ASSERT(!isnan(tt_resolution(buffer)),NULL);
    
    //Make sure that there is data - and if not, return empty result
    TT_WARN(tt_datanum(buffer)==0,return (coincidenceMatrix?coincidenceMatrix:(uint64_t*)calloc(tt_channels(buffer)*tt_channels(buffer),sizeof(uint64_t)));,"Buffer is empty!");

    //Convert the delay array to bins
    if (delayArray) {
        delays = tt_delaytime2bin(buffer,delayArray,NULL,tt_channels(buffer));
        TT_WARN(!delays,,"Delays array allocation failed. Ignoring delays.");
    }
    
    coincidenceMatrix = tt_rawcoincidences(buffer, tt_subtractreference(buffer,tt_time2bin(buffer,time)), tt_time2bin(buffer,radius), 
            coincidenceMatrix, delays, tt_datanum(buffer)-1);
    
    //Free the allocated delay array
    tt_free(delays);
    
    return coincidenceMatrix;
}

TT_DEF_ uint64_t tt_rawmulticoincidences(const tt_buf *const buffer, uint64_t timebins, uint64_t diameter, uint8_t* channels, uint8_t channelnum,
        uint64_t* delayArray, uint64_t dataindex) {
    uint64_t* channelIndex; 
    uint64_t* channelMin;
    uint64_t* channelMax;
    uint64_t icounter;
    uint64_t coincidences = 0;
    int channelCoincidence;
    int i,j;
    
    //Make sure everything is working correctly
    //TT_ASSERT(buffer,NULL);
    TT_ASSERT(timebins,0);
    TT_ASSERT(channels,0);
    TT_ASSERT(channelnum,0);
    TT_ASSERT(dataindex < tt_datanum(buffer),0);
    
    //Allocate all necessary memory
    
    //This all can be turned into one memory allocation, at the cost of readability
    channelIndex = (uint64_t*)malloc(sizeof(uint64_t)*channelnum);
    TT_CHKERR(channelIndex,return 0,"Memory allocation failed!");
    channelMin = (uint64_t*)malloc(sizeof(uint64_t)*channelnum);
    TT_CHKERR(channelMin,free(channelIndex);return 0,"Memory allocation failed!");
    channelMax = (uint64_t*)malloc(sizeof(uint64_t)*channelnum);
    TT_CHKERR(channelMax,free(channelMin); free(channelIndex); return 0,"Memory allocation failed!");
    
    //Find channel iterator locations, and set the bounds for channel times (channelMin, channelMax)
    for (i=0; i< channelnum;i++) {
    
        //put the index at current location
        channelIndex[i] = dataindex;
        
        //Set the channel's maximum timing
        channelMax[i] = tt_tag(buffer,dataindex);
        
        //If there is a delay array, find iterator locations offset by the delay amounts
        if (delayArray) {
            channelIndex[i] -= tt_bins2points(buffer,dataindex,delayArray[i]);
            //printf("InitialIndex: %i %i %i %i %i\n",i,(int)channels[i],(int)channelIndex[i],(int)tt_tag(buffer,channelIndex[i]),(int)tt_channel(buffer,channelIndex[i]));
            //The difference between bins2points and the intended result is that they behave differently when there are multiple points with the same time tag
            //  on the boundary
            //For example, with boundary being 1200,  the datapoint returned is:
            //         100,1000,1200,1200,1200,1200,1400,1800
            //resut:            ^
            //want:                            ^
            //These are, of course, details, but they are important nevertheless! We have to account for them!
            for (icounter=channelIndex[i];icounter < dataindex && tt_tag(buffer,icounter+1)==tt_tag(buffer,channelIndex[i]);icounter++);
            channelIndex[i]=icounter;
            
            //printf("InitialIndex^: %i %i %i %i %i\n",i,(int)channels[i],(int)channelIndex[i],(int)tt_tag(buffer,channelIndex[i]),(int)tt_channel(buffer,channelIndex[i]));
            //Make sure that the channel's Max time does not overflow. Look up TT_WARN to find
            //  out why putting an else there is kosher.
            TT_WARN(channelMax[i] <= delayArray[i],channelMax[i]=0,"Delay for channel %i goes beyond available time!",(int)channels[i]) 
                else {
                    channelMax[i] -= delayArray[i];
                }
        }
        
        //Next, subtract the allowed time to get the minimum time tag to accept. The minindex is taken into account later
        TT_WARN(channelMax[i] <= timebins, channelMin[i] = 0,"Time for channel %i goes beyond total data!",(int)channels[i])
            else {
                channelMin[i] = channelMax[i] - timebins;
            }
    }
    
    //Now make sure that the channelIndices are pointing to the correct channels, and that they are within bounds
    for (i=0;i< channelnum;i++) {
    
        //Make sure that the channel index is within channelMax
        while (((~((uint64_t)0))!=channelIndex[i]) && tt_tag(buffer,channelIndex[i]) > channelMax[i]) {
            channelIndex[i] = ttinternal_getNextChannelIndex(buffer,channelIndex[i], channels[i],channelMin[i],tt_minindex(buffer));
        }
        
        
        if (((~((uint64_t)0))!=channelIndex[i]) && channels[i]!=(int)tt_channel(buffer,channelIndex[i])) {
            //Not the correct channel, move index pointer to the correct one.
            
            channelIndex[i] = ttinternal_getNextChannelIndex(buffer,channelIndex[i], channels[i],channelMin[i],tt_minindex(buffer));
        }
        
        //Make sure each channel has an index. If not, then there is no use in finding coincidences
        if (~((uint64_t)0)==channelIndex[i]) {
            
            //Free all necessary memory
            free(channelIndex);
            free(channelMin);
            free(channelMax);
            return 0;
        }
    }
    //for (i=0;i<channelnum;i++) {
    //    printf("DELAY: %i %i MAX:%i\n",(int)channels[i],(int)delayArray[i],(int)channelMax[i]);
    //}
    //for (i=0;i<channelnum;i++) {
    //    printf("INDEX %i %i %i\n",(int)channels[i],(int)channelIndex[i],(int)tt_tag(buffer,channelIndex[i]));
    //}
    //Preliminaries are set. This means:
    // - channelIndex contains the index of the first element of each channel
    // - channelMax contains the time from which to count down 'timebins' for each time
    // - channelMin is the minimum time tag to allow for each channel
    
    while ((i= ttinternal_getMinChannel(buffer,channelMax,channelIndex,channelnum))!=-1) {
        //Set the number of channels in coincidence to 0
        channelCoincidence = 0;
        //printf("CMP: %i %i %i %i\n",i,channels[i],(int)channelIndex[i],(int)tt_tag(buffer,channelIndex[i]));
        //Take the smallest, and check if it has a coincidence with ALL others
        for (j=0; j< i; j++) {            
            //printf("F %i %i | %i %i | %i %i\n",(int)channels[i],(int)channels[j],(int)channelIndex[i],(int)channelIndex[j],(int)tt_tag(buffer,channelIndex[i]),(int)tt_tag(buffer,channelIndex[j]));
            
            if (!(channelMax[i]-tt_tag(buffer,channelIndex[i])+diameter >= channelMax[j]-tt_tag(buffer,channelIndex[j]))) {
                j=-1;   //Set j to -1, so that we know that it is invalid
                break;
            }
            else {
                //printf("Channel %i coincidence\n",channels[j]);
                channelCoincidence++;
            }
        }
        if (j!=-1) {
            for (j=i+1; j< channelnum; j++) {            
                //printf("L %i %i | %i %i | %i %i\n",(int)channels[i],(int)channels[j],(int)channelIndex[i],(int)channelIndex[j],(int)tt_tag(buffer,channelIndex[i]),(int)tt_tag(buffer,channelIndex[j]));
                
                if (!(channelMax[i]-tt_tag(buffer,channelIndex[i])+diameter >= channelMax[j]-tt_tag(buffer,channelIndex[j]))) {
                    break;
                }
                else {
                    //printf("Channel %i coincidence\n",channels[j]);
                    channelCoincidence++;
                }
            }
        }
        if (channelCoincidence==channelnum-1) coincidences++;
        
        channelIndex[i] = ttinternal_getNextChannelIndex(buffer,channelIndex[i], channels[i], channelMin[i],tt_minindex(buffer));
        
        //If any of the channels goes out of bounds, then by definition, we cannot get another coincidence
        if (~((uint64_t)0)==channelIndex[i]) break;
    }
    
    //Free all necessary memory
    free(channelIndex);
    free(channelMin);
    free(channelMax);
    
    return coincidences;
}

//Multiple concidences - says that there was a coincidence only if all channels of array present within diameter
TT_DEF_ uint64_t tt_multicoincidences(const tt_buf *const buffer, double time, double diameter, uint8_t* channels, int channelnum, double* delayArray) {
    uint64_t* delays = NULL;
    uint64_t result;
    
    //TT_ASSERT(buffer,0);
    TT_ASSERT(time,0);
    
    //Check for reference time
    TT_ASSERT(!isnan(tt_resolution(buffer)),0);
    
    //Make sure that there is data - and if not, return empty result
    TT_WARN(tt_datanum(buffer)==0,return 0;,"Buffer is empty!");

    if (delayArray) {
        delays = tt_delaytime2bin(buffer,delayArray,NULL,channelnum);
        TT_WARN(!delays,,"Delays array allocation failed. Ignoring delays.");
    }
    
    result =  tt_rawmulticoincidences(buffer, tt_subtractreference(buffer,tt_time2bin(buffer,time)), 
            tt_time2bin(buffer,diameter),channels, (uint8_t)channelnum, delays,tt_datanum(buffer)-1);
    
    tt_free(delays);
    
    return result;
}

TT_DEF_ uint64_t* tt_rawcorrelate(const tt_buf *const buffer, uint64_t timebins, uint64_t windowradius, int bins, uint8_t channel1, uint64_t delay1, \
                                    uint8_t channel2, uint64_t delay2, uint64_t* resultArray, uint64_t dataindex) {
    //This has two iterators, one chasing the other. The first goes one way and just finds each element in the first channel, and the second follows, going into
    //  the windowtime of the first, and adding itself to coincidences. After iterating once, it moves the iterator back to the edge of windowtime for the second sequence,
    //  and does the same thing
    
    //For speed reasons, can later add a buffer for the second channel of ~20 points, so that moving back is faster.
    
    uint64_t channelMin[2];
    uint64_t channelIndex[2];
    uint64_t channelMax[2];
    uint64_t maxIndex;
    uint64_t icounter,lowerBound,upperBound;
    uint64_t curindex;
    double clicksperbin;
    //TT_ASSERT(buffer,NULL);
    TT_ASSERT(timebins,NULL);
    TT_ASSERT(windowradius,NULL);
    TT_ASSERT(bins, NULL);
    TT_ASSERT(dataindex < tt_datanum(buffer),NULL);    //If this is wrong, can cause infinite looping
    
    clicksperbin = (windowradius*2)/(double)bins;
    TT_ASSERT(clicksperbin>=1.0,NULL);    //Make sure that each bin is at least one time tag
    
    //Now, do the same thing as in both coincidence functions to find the bounds and indices for the two channels
    //Find channel iterator locations, and set the bounds for channel times (channelMin, channelMax)
    
    channelIndex[0] = dataindex;                                            //put the index at current location
    channelMax[0] = tt_tag(buffer,dataindex);                               //Set the channel's maximum timing
    channelIndex[0] -= tt_bins2points(buffer,dataindex,delay1);             //Correct for delays
    channelIndex[1] = dataindex;                                            //put the index at current location
    channelMax[1] = tt_tag(buffer,dataindex);                               //Set the channel's maximum timing
    channelIndex[1] -= tt_bins2points(buffer,dataindex,delay2);             //Correct for delays
    
    //WARNING: This should never happen
    TT_CHKERR(channelIndex[0] <= dataindex && channelIndex[1] <= dataindex,return NULL,"Impossibility happened. Is the input data ordered?");
    
    //First, create the result array if such was not yet created
    if (!resultArray) {
        resultArray = (uint64_t*)calloc(sizeof(uint64_t),bins);
        TT_CHKERR(resultArray,return NULL,"Memory allocation failed!");
    }
    
    //The difference between bins2points and the intended result is that they behave differently when there are multiple points with the same time tag
    //  on the boundary
    //For example, with boundary being 1200,  the datapoint returned is:
    //         100,1000,1200,1200,1200,1200,1400,1800
    //resut:            ^
    //want:                            ^
    //These are, of course, details, but they are important nevertheless! We have to account for them!
    for (icounter=channelIndex[0];icounter < dataindex && tt_tag(buffer,icounter+1)==tt_tag(buffer,channelIndex[0]);icounter++);
    channelIndex[0]=icounter;
    for (icounter=channelIndex[1];icounter < dataindex && tt_tag(buffer,icounter+1)==tt_tag(buffer,channelIndex[1]);icounter++);
    channelIndex[1]=icounter;

    
    //Make sure that the channel's Max time does not overflow. Look up TT_WARN to find
    //  out why putting an else there is kosher.
    TT_WARN(channelMax[0] <= delay1,channelMax[0]=0,"Delay for channel %i goes beyond available time!",(int)channel1) 
        else channelMax[0] -= delay1; 
    TT_WARN(channelMax[1] <= delay2,channelMax[1]=0,"Delay for channel %i goes beyond available time!",(int)channel2) 
        else channelMax[1] -= delay2;
    
    //Next, subtract the allowed time to get the minimum time tag to accept. The minindex is taken into account later
    TT_WARN(channelMax[0] <= timebins, channelMin[0] = 0,"Time for channel %i goes beyond total data!",(int)channel1)
        else channelMin[0] = channelMax[0] - timebins;
    TT_WARN(channelMax[1] <= timebins, channelMin[1] = 0,"Time for channel %i goes beyond total data!",(int)channel2)
        else channelMin[1] = channelMax[1] - timebins;
    
    //Put the channelIndex within bounds.
    while (((~((uint64_t)0))!=channelIndex[1]) && tt_tag(buffer,channelIndex[1]) > channelMax[1]) {
        channelIndex[1] = ttinternal_getNextChannelIndex(buffer,channelIndex[1], channel2,channelMin[1],tt_minindex(buffer));
    }
    
    if (((~((uint64_t)0))!=channelIndex[1]) && channel2!=(int)tt_channel(buffer,channelIndex[1])) {
        channelIndex[1] = ttinternal_getNextChannelIndex(buffer,channelIndex[1], channel2,channelMin[1],tt_minindex(buffer));
    }
    
    //Make sure each channel has an index. If not, then the correlation is all 0s,
    if (~((uint64_t)0)==channelIndex[1]) return resultArray;
    
    //Put the channelIndex within bounds.
    while (((~((uint64_t)0))!=channelIndex[0]) && tt_tag(buffer,channelIndex[0]) > channelMax[0]) {
        channelIndex[0] = ttinternal_getNextChannelIndex(buffer,channelIndex[0], channel1,channelMin[0],tt_minindex(buffer));
    }
    
    //Now make sure that the channelIndices are pointing to the correct channels
    
    if (((~((uint64_t)0))!=channelIndex[0]) && channel1!=(int)tt_channel(buffer,channelIndex[0])) {
        //Not the correct channel, move index pointer to the correct one.
        channelIndex[0] = ttinternal_getNextChannelIndex(buffer,channelIndex[0], channel1,channelMin[0],tt_minindex(buffer));
    }
    
    
    //We don't want channel 2 to go above the max index for channel 2 ever. This only applies to channel 2.
    maxIndex = channelIndex[1];
    
    //Make bins one smaller, since it is used as an index marker of the last element later
    --bins;
    
    //The bounds are set for both channels. Now it is time to iterate through the elements of the first channel, one by one.
    while (~((uint64_t)0)!=channelIndex[0]) {
        
        //First, find the bounds of the "VALID REGION" - the region around channelIndex[0] where channelIndex[1] can be put.
        //  The conversions are to put it in terms of channel 2's time tags.
        
        //upperBound = channelMax[1] - (channelMax[0]-(tt_tag()+windowradius))
        //upperBound = channelMax[1]+tt_tag(buffer,channelIndex[0])+windowradius-channelMax[0];
        upperBound = (channelMax[1]+tt_tag(buffer,channelIndex[0])+windowradius > channelMax[0] ? channelMax[1]+tt_tag(buffer,channelIndex[0])+windowradius-channelMax[0] : 0);
        //lowerBound = channelMax[1] - (channelMax[0]-(tt_tag()-windowradius))
        //lowerBound = channelMax[1]+tt_tag(buffer,channelIndex[0])-windowradius-channelMax[0];
        lowerBound = (channelMax[1]+tt_tag(buffer,channelIndex[0]) > channelMax[0]+windowradius ? channelMax[1]+tt_tag(buffer,channelIndex[0])-windowradius-channelMax[0] : 0);
        
        //If the channel's minimum is greater than the upper bound, exit loop
        if (upperBound <= channelMin[1]) break;
        //If the minimum is greater than the lower bound, set lowerBound to minimum
        if (lowerBound < channelMin[1]) lowerBound = channelMin[1];
        
        //Move the second channel cursor back to beyond the upper bound if it is within bounds
        while (maxIndex > channelIndex[1] && tt_tag(buffer,channelIndex[1]) <= upperBound) {
             ++channelIndex[1];
        }

        //Move the second channel cursor forward to the edge of the window size (not necessarily the correct channel), meaning that
        //  this datapoint is right INSIDE the bounds
        while (tt_minindex(buffer) < channelIndex[1] && tt_tag(buffer,channelIndex[1]) > channelMin[1] && tt_tag(buffer,channelIndex[1]) > upperBound) {
            --channelIndex[1];
        }
        
        if ((tt_tag(buffer,channelIndex[1]) <= upperBound)) {
            //At this point in time, channelIndex[1] points to the first data point within upper bound - though not necessarily of the correct channel
            
            //Check if the current data point is of the correct channel, and within bounds. If it is, add it to the result at the correct location
            if (tt_channel(buffer,channelIndex[1]) == channel2 && tt_tag(buffer,channelIndex[1]) > lowerBound) {
                TT_WARN((upperBound-tt_tag(buffer,channelIndex[1]))/clicksperbin > bins,,"+INDEX FAILURE: %f\n",((double)(upperBound-tt_tag(buffer,channelIndex[1])))/clicksperbin)
                else
                    ++resultArray[bins-(int)((upperBound-tt_tag(buffer,channelIndex[1]))/clicksperbin)];
            }

            //WARNING:There might be a negative index here!!
            while ((curindex = ttinternal_getNextChannelIndex(buffer,channelIndex[1],channel2,lowerBound+1,tt_minindex(buffer)))!=~((uint64_t)0)) {
                channelIndex[1] = curindex;
                //The given tag is within bounds of the binning, so find which bin it belongs in, and increment it!
                //Since the lower bound is not guaranteed to exist, we use the upper bound as a stopping point
                TT_WARN(((int)((upperBound-tt_tag(buffer,channelIndex[1]))/clicksperbin)>bins),,"INDEX FAILURE: %f\n",(upperBound-tt_tag(buffer,channelIndex[1]))/clicksperbin)
                else
                    ++resultArray[bins-(int)((upperBound-tt_tag(buffer,channelIndex[1]))/clicksperbin)];
            }
        }
        //Get the next index of the first channel
        channelIndex[0] = ttinternal_getNextChannelIndex(buffer,channelIndex[0], channel1, channelMin[0],tt_minindex(buffer));
    }

    return resultArray;
}

//time - the amount of data to go through
//windowradius - the time window to correlate
//bins - the number of bins +- to divide the window. The resulting array will have bins elements
TT_DEF_ uint64_t* tt_correlate(const tt_buf *const buffer, double time, double windowradius, int bins, uint8_t channel1, double delay1, uint8_t channel2,  double delay2,uint64_t* resultArray) {
    double delayArray[2];
    uint64_t delays[2];
    
    //TT_ASSERT(buffer,0);
    TT_ASSERT(time,0);
    
    //Check for reference time
    TT_ASSERT(!isnan(tt_resolution(buffer)),NULL);
    
    //Make sure that there is data - and if not, return empty result
    TT_WARN(tt_datanum(buffer)==0,return (resultArray?resultArray:(uint64_t*)calloc(sizeof(uint64_t),bins));,"Buffer is empty!");

    //Create delays
    delayArray[0] = delay1;
    delayArray[1] = delay2;
    //The function should always return the array passed in
    TT_ASSERT(delays==tt_delaytime2bin(buffer,delayArray,delays,2),NULL);
    
    //Debug printout: subtractreference ws giving weird results unpredictably
    //printf("CORRELATE:\nTIME: %f\nwrad: %e\nBINS: %i\nD: %e %e\n",time,windowradius,bins,delay1,delay2);
    //printf("\nD1: %llu\nD2: %llu\nTIME: %llu\n",delays[0],delays[1],tt_subtractreference(buffer,tt_time2bin(buffer,time)));
    return tt_rawcorrelate(buffer,tt_subtractreference(buffer,tt_time2bin(buffer,time)),
            tt_time2bin(buffer,windowradius),bins,channel1,delays[0],channel2,delays[1],resultArray,tt_datanum(buffer)-1);
}
