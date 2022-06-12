'''
By Andreas H. Jørgensen
'''
'''
The goal of this file is to verify that the probability that the
k=m/3 largest samples out of m sample from the rounded gaussian/normal
distribution is at p/2trk with high (non-negligible probability).
This is to verify that for any row in the error matrix E of the cryptosystem,
we will have that the k=m/3 largest elements sums to at least p/2tr sucht that
when scaling each of the k elements by r we have a sum of at least p/2t giving
rise to a decryption error.
-> this gives os the possibility to chose d such that there is a high probability
   of decryption errors

'''

import numpy as np
import time
import multiprocessing as mp

def sample(m,p,alpha):
    '''
    Used to make m samples from discrete Psi_alpha distribution,
    with standard deviation p*alpha/sqrt(2pi)
    '''
    std = alpha*p/np.sqrt(2*np.pi)
    errors = np.mod(np.random.normal(0,std,m).round(),p)
    return errors

def canGiveDecErr(errs,p,t,r):
    '''
    Takes m error-values and determines if all of the largest k=ceil(m/3) absolute values are
    at least p/2trk.
    '''
    m,=errs.shape
    k=int(np.ceil(m/3))
    errs[errs>p/2]=p-errs[errs>p/2]
    sorted = np.sort(errs)[::-1]
    return sorted[k-1]>=p/(2*t*r*k)

def testError(m,p,t,r,alpha):
    '''
    Uses the set of parameters to generate m samples (a row of E)
    and checks whether we can use it to enforce an encryption error
    Using LWE oracle
    '''
    errs = sample(m,p,alpha)
    res = canGiveDecErr(errs,p,t,r)
    return res
    # if(canGiveDecErr(errs,p,t,r)):
        # return 1
    # else:
        # return 0
    #

def errorProb(m,p,t,r,alpha,repetitions,elaborate=False):
    '''
    Uses the set of parameters to generate m samples
    and checks whether we can use it to enforce an encryption error
    Using LWE oracle.
    This is repeated 'repetitions' times and a success probability is
    estimated and returned
    '''
    count = 0
    for i in range(repetitions):
        count += testError(m,p,t,r,alpha)
        if((i+1)%2000==0):
            print("(" + str(i+1)+ "/" + str(repetitions)+ "=" + "{:.1f}".format(i/repetitions*100)+"%) Success Rate: " + "{:.2f}".format(count/(i+1)*100) + "%")
    return count/repetitions

def errorProbParallel(m,p,t,r,alpha,repetitions,elaborate=False):
    '''
    Parallelized version of errorProb
    Uses the set of parameters to generate m samples
    and checks whether we can use it to enforce an encryption error
    Using LWE oracle.
    This is repeated 'repetitions' times and a success probability is
    estimated and returned
    '''
    
    with mp.Pool(processes=mp.cpu_count()+1) as pool:
        results = list()
        for i in range(repetitions):
            res = pool.apply_async(testError,(m,p,t,r,alpha))
            results.append(res)
        
        count = 0        
        for i,res in enumerate(results):
            #print(res.get())
            #if res.get()==True:
            #    count = count+1
            count += res.get()
            if((i+1)%2000==0):
                print("(" + str(i+1)+ "/" + str(repetitions)+ "=" + "{:.1f}".format(i/repetitions*100)+"%) Success Rate: " + "{:.2f}".format(count/(i+1)*100) + "%")
    return count/repetitions*100

def main():
    #m=10;p=10;t=2;r=5;alpha=0.65
    m=4536;p=32749;t=40;r=1;alpha=0.000217
    #m=2008;p=2003;t=2;r=1;alpha=0.0065
    #m=1500;p=8191;t=4;r=5;alpha=0.0009959
    #m=1319;p=4093;t=2;r=4;alpha=0.0024
    
    #My own parameters (See word)
    #n=136; l=136; m=8000;  p=34000; r=1; t=2; alpha=0.006321
    #n=136; l=136; m=8001;  p=34000; r=1; t=2; alpha=0.0053
    #n=136; l=136; m=8002;  p=34000; r=1; t=2; alpha=0.001528265
    
    totalCount = 100000
    count = 0
    for i in range(totalCount):
        errs = sample(m,p,alpha)
        res = canGiveDecErr(errs,p,t,r)
        #print(res)
        if(res==True):
            count+=1
        if(i%2000==0):
            print("(" + str(i)+ "/" + str(totalCount)+ "=" + "{:.1f}".format(i/totalCount*100)+"%) Success Rate: " + "{:.2f}".format(count/(i+1)*100) + "%")
    print("(Complete) Success Rate:" + str(count/totalCount*100) + "%")

if __name__=="__main__":
    main()