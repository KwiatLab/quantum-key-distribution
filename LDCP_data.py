
from numpy import *
from SW_prep import *
from entropy_calculator import *
import sys
alice=load("resultsLaurynas/aliceChannels_06032014_maxpower.npy")
bob=load("resultsLaurynas/bobChannels_06032014_maxpower.npy")

#Assumes the person's data is ordered
def createLDPCdata(person,tot=None,alphabet=16):
    
    sw_p=person.copy()
    sw_p/=alphabet
    
    if (tot==None): tot=sw_p[-1]+1
    
    person_sw=zeros(tot,dtype=int8)
    person_nb = zeros(tot,dtype=int32)
    
    doublenum=0
    for z in xrange(len(person)):
        if (person_sw[sw_p[z]]!=0):
            if (person_sw[sw_p[z]]==1):
                doublenum+=1
            #There is already a value here, delete it - this is the 2-pulse rate
            person_sw[sw_p[z]]+=1
        else:
            #There is no value here, so save the value
            person_sw[sw_p[z]]=1
            person_nb[sw_p[z]]=person[z]%alphabet
    
    return (person_sw,person_nb,doublenum)
def probLetter(l,alph):
    p=zeros(alph)
    for i in xrange(alph):
        p[i]=sum(l==i)
    p/=len(l)
    return p
    

for alphabet in [128]:#xrange(201,1000,50): #array(range(1,10)):
    print "RUNNING",alphabet
    tot = alice[-1]/alphabet
    if (tot< bob[-1]/alphabet): tot = bob[-1]/alphabet
    tot+=1
    print "->Creating LDPC arrays"
    (alice_sw,alice_nb,alice_doublenum) = createLDPCdata(alice,tot,alphabet)
    (bob_sw,bob_nb,bob_doublenum)=createLDPCdata(bob,tot,alphabet)
    print "\tDoubles (Alice):",alice_doublenum,float(alice_doublenum)/len(alice_sw)
    print "\tDoubles (Bob):",bob_doublenum,float(bob_doublenum)/len(bob_sw)
    
    print "->Extracting ideal NB-LDPC arrays"
    sys.stdout.flush()
    
    nb_mask = logical_and(alice_sw,bob_sw)
    
    
    alice_nbf = alice_nb[nb_mask]
    bob_nbf = bob_nb[nb_mask]
    
    print "SAVING FILE"
    savetxt("FRAME_128_DATA.csv",(alice_nbf,bob_nbf),fmt="%i")
    
    #The probability of a 1 in the original-original sequence
    p1 = float(len(alice))/alice[-1]
    print "p1:",p1,float(len(bob))/bob[-1]
    
    total_c = intersect1d(alice,bob)
    p1g1 = float(len(total_c))/len(alice)
    print "Coincidence rate (p1g1):",p1g1,float(len(total_c))/len(bob)
    
    if (any(alice_sw >= alphabet) or any(bob_sw >= alphabet)):
        print "WARNING: Over the TOP!"
        alice_sw[alice_sw >= alphabet]=alphabet-1
        bob_sw[bob_sw >= alphabet]=alphabet-1
    swtransmat = transitionMatrix_data2(alice_sw,bob_sw,alphabet)
#if (any(transitionMatrix_data2(alice_sw,bob_sw,alphabet)!=swtransmat)):
    #    print "NOT THE SAME!!!!!\n\n\n"
    swpl = probLetter(alice_sw,alphabet)
    print "Letter Probabilities:"
    print swpl
    print "Transition Matrix (SW):"
    print swtransmat
    
    nbtransmat = transitionMatrix_data2(alice_nbf,bob_nbf,alphabet)
    nbpl = probLetter(alice_nbf,alphabet)
    print "Letter Probabilities:"
    print nbpl
    print "Transition Matrix (NB):"
    print nbtransmat
    
    #The total number of original bins is alice, and the final bits is the number of nonbinary left
    nb_bperf= alice[-1]/len(bob_nbf)
    print "Number of original bits per nonbinary:",nb_bperf
    
    (te,te2,be,nbe)=entropy_calculate2(p1,p1g1,p1,0.34,alphabet,swpl,swtransmat,nb_bperf,nbpl,nbtransmat)
    f=open("RESULT_2013_11_2.csv","a")
    f.write(str(alphabet)+" "+str(te)+" "+str(te2)+" "+str(be)+" "+str(nbe)+"\n")
    f.close()
    
