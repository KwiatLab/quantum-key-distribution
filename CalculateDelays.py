'''
Created on Feb 28, 2016

@author: laurynas
'''

from pylab import *
#Fit to the given curve
from numpy import *
from scipy.optimize import curve_fit
import ttag
import graphs

from itertools import product

rc('text', usetex=True)
rc('font', family='serif')

# Define model function to be used to fit to the data above:

def gauss(x, *p):
    A, mu, sigma = p
    return A*exp(-(x-mu)**2/(2.*sigma**2))


#Assume that the delay is < 1ms
#Also, just assume that max one is the correct delay - a gaussian fit might be better
#   for certain things
def getDelay(bufAlice,channel1,channel2,initialdelay1=0.0,initialdelay2=0.0,delaymax = 0.0000001,time=1.0):
    bins = int(delaymax/bufAlice.resolution)*2
    corr = bufAlice.correlate(time,delaymax,bins,channel1,channel2,channel1delay=initialdelay1,channel2delay=initialdelay2)
    
    #Now, we have a way to fit to gaussian - set initial parameters
    mu = argmax(corr)
    sigma = 5
    A = max(corr)
    popt,pcov = curve_fit(gauss,range(bins),corr,p0=(A,mu,sigma))
    print(channel1,channel2,"FIT: (A,mu,sigma)=",popt)
    return (popt[1]-len(corr)/2)*bufAlice.resolution


#This function cannot be used on huge buffers, since it creates a copy of the entire dataset
def getPossibleInitialDelays(bufAlice,syncChannel1,syncChannel2):
    channels,singles = bufAlice[:]
    
    arr1 = where(channels==syncChannel1)
    arr2 = where(channels==syncChannel2)
    print arr1
    print arr2
    possibilities1= arr1[0]
    possibilities2= arr2[0]

    delays=[]

    for i1,i2 in product(possibilities1,possibilities2):
        delays.append(singles[i2]-singles[i1])

    return delays

def getDelays(bufAlice,channels1,channels2,initialdelay2=0.0,delays1=None,delays2=None,delaymax=0.0000001,time=1.0):

    if (delays1==None):
        delays1=zeros(len(channels1))
    if (delays2==None):
        delays2=ones(len(channels2))*initialdelay2

    #First set all of channels2 delays
    for i in range(len(delays2)):
        delays2[i] += getDelay(bufAlice,channels1[0],channels2[i],delays1[0],delays2[i],delaymax=delaymax,time=time)

    #Next, set all of delays for channels1
    for i in range(1,len(delays1)):
        delays1[i] -= getDelay(bufAlice,channels1[i],channels2[0],delays1[i],delays2[0])

    return (delays1,delays2)
"""
def plotAll(bufAlice,channels1,channels2,delays1,delays2):
    bins = 100
    length = bins/2*bufAlice.resolution
    dist=(array(range(bins))-bins/2)*bufAlice.resolution
    f,ax=subplots(len(channels1),sharex=True)
    plots = []
    for i in range(len(channels1)):
        for j in range(len(channels2)):
            ax[i].set_ylabel("Channel "+str(i+1))
            cor=bufAlice.correlate(1.0,length,bins,channels1[i],channels2[j],channel1delay=delays1[i],channel2delay=delays2[j])
            #Now, we have a way to fit to gaussian - set initial parameters
            mu = argmax(cor)
            sigma = 5
            A = max(cor)
            popt,pcov = curve_fit(gauss,range(bins),cor,p0=(A,mu,sigma))
            plots.append(ax[i].plot(dist,cor,linewidth=2,label=r"Channel "+str(j+1) + r" (A="+str(int(round(popt[0])))+r" $\sigma$="+str(round(popt[2]*bufAlice.resolution*1e12,2))+"ps)"))
    ax[0].set_title("Correlations Between Alice and Bob's Channels")
    f.subplots_adjust(hspace=0)
    ax[0].legend()
"""

