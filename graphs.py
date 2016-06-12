from pylab import *
#Fit to the given curve
from numpy import *
from scipy.optimize import curve_fit
import ttag

def gauss(x, *p):
    A, mu, sigma = p
    return A*exp(-(x-mu)**2/(2.*sigma**2))

def plotAutocorrelations(bufAlice,time=1.0,bins=2000,ymax=None):
    length=bins/2*bufAlice.resolution
    dist=(array(range(bins))-bins/2)*bufAlice.resolution
    f,ax=subplots(bufAlice.channels,sharex=True)
    for i in range(bufAlice.channels):
        ax[i].set_ylabel("Channel "+str(i+1))
        if (ymax):
            ax[i].set_ylim((0,ymax))
            ax[i].set_xlim((dist[0],dist[-1]))
        ax[i].plot(dist,bufAlice.correlate(time,length,bins,i,i))
    ax[0].set_title("Autocorrelations Of Each Channel")
    ax[-1].set_xlabel("Time in Seconds")
    f.subplots_adjust(hspace=0)
    show()

def plotABCorrelations(bufAlice,channels1,channels2,delays1=None,delays2=None,bins=150,pulsebin=0.0,time=1.0):
    if (delays1==None):
        delays1=zeros(len(channels1))
    if (delays2==None):
        delays2=zeros(len(channels2))
        
    window_radius = bins/2*bufAlice.resolution
    dist=(array(range(bins))-bins/2)*bufAlice.resolution
    f,ax=subplots(len(channels1),sharex=True)
    plots = []
    for i in range(len(channels1)):
        ax[i].set_ylabel("Channel "+str(channels1[i]))
        if (pulsebin > 0.0):
            ax[i].axvline(pulsebin/2)
            ax[i].axvline(-pulsebin/2)
            ax[i].axvline(3*pulsebin/2)
            ax[i].axvline(-3*pulsebin/2)
        for j in range(len(channels2)):
            cor=bufAlice.correlate(time,window_radius,bins,channels1[i],channels2[j])
            print "This is correlations between "+str(channels1[i])+" and "+str(channels2[j])+" : \n",cor, max(cor)
#             #Now, we have a way to fit to gaussian - set initial parameters
            popt=[max(cor),argmax(cor),3]
            try:
                mu = argmax(cor)
                sigma = 5
                A = max(cor)
                popt,pcov = curve_fit(gauss,range(bins),cor,p0=(A,mu,sigma))
            except:
                popt=[max(cor),argmax(cor),3]
#             #
            if (i==0):
                plots.append(ax[i].plot(dist,cor,linewidth=2,label=r"Channel "+str(j+1) + r" ($\sigma$="+str(round(popt[2]*bufAlice.resolution*1e12,2))+"ps)"))
            else:
                plots.append(ax[i].plot(dist,cor,linewidth=2,label=r"$\sigma$="+str(round(popt[2]*bufAlice.resolution*1e12,2))+"ps"))
            ax[i].legend()
    ax[0].set_title("Correlations Between Alice and Bob's Channels")
    ax[-1].set_xlabel("Time in Seconds")
    f.subplots_adjust(hspace=0)
    ax[0].legend()
    show()

def polarizationPlot(polA,polB):
    figure()
    suptitle("Polarization Channels")
    for i in range(4):
        subplot(2,2,i+1)
        title("Channel "+str(1+i))

        v = polA[polB==i]

        bar([1,2,3,4],array([sum(v==0),sum(v==1),sum(v==2),sum(v==3)],dtype=double)/float(len(v)))
        #if (i==0):
        #    legend()
    show()


if (__name__=="__main__"):
    aliceChannels=[0,1,2,3]
    bobChannels=[4,5,6,7]
    binsize = 1/(32*120.0e6)
    bufAlice = ttag.TTBuffer(1)
    plotAutocorrelations(bufAlice,ymax=3000,bins=10000)
    plotABCorrelations(bufAlice,aliceChannels,bobChannels,bins=60,pulsebin=binsize)