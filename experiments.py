from LWE import LWE
import numpy as np
from scipy.stats import norm

elaborate=True
n=136; l=136; m=8000;  p=34000; r=1; t=2; alpha=0.006321
#                   n       l       m       p       r       t   alpha       epsilon     error Prob,
parameters = [  (   136,    136,    8000,   34000,  1,      2,  -1,         0.01,       -1), #neglible error
                (   136,    136,    8001,   34000,  1,      2,  -1,         0.1,        -1), #negliblible error
                (   136,    136,    8002,   34000,  1,      2,  -1,         0.1,        0.01), #1% error
                (   136,    136,    2008,   2003,   1,      2,  0.0065,     -1,         -1), #Examples from book chapter (does not satisfy p>=2tmr)
                (   166,    166,    1319,   4093,   4,      2,  0.0024,     -1,         -1),
                (   192,    192,    1500,   8191,   5,      4,  0.0009959,  -1,         -1),
                (   214,    214,    1333,   16381,  12,     4,  0.00045,    -1,         -1),
                (   233,    233,    1042,   32749,  59,     2,  0.000217,   -1,         -1),
                (   233,    233,    4536,   32749,  1,      40, 0.000217,   -1,         -1)
             ]

parameters = []             

for i,(n,l,m,p,r,t,alpha, epsilon, pr) in enumerate(parameters):
    if alpha==-1:
        if pr==-1:
            if epsilon == -1:
                print("ERROR: both, alpha, epsilon and error probility is undifined for", i,"set of parameter")
                continue #Go to next set of parameters
            g=np.log2(n)**(1/2+epsilon) #IMPOTANT: MAYBE CHANGE TO LOG INSTEAD OF LOG2
            #We find alpha based on goal, that we want a negligible probability of decryption errors
            alpha = 1/(t*np.sqrt(r*(r+1)*m/(6*np.pi))*g)
            #We also find the probability of decryption errors
            pr = 4/g*np.exp(-(g/4)**2/2) #IMPORTANT: SHOULD BE UPDATED BASED ON NEW BOUND ON ROUNDING ERROR
        else:
            #We are given a specific target probability of decryption errors, and we find alpha based on this
            quantile = norm.ppf(1-pr/2)
            alpha = 1/(4*t*np.sqrt(r*(r+1)*m/(6*np.pi))*quantile) #IMPORTANT: SHOULD BE UPDATED BASED ON NEW BOUND ON ROUNDING ERROR
    else:
        #We already know alpha
        #We will now find epsilon
        epsilon=np.log(1/(t*np.sqrt(r*(r+1)*m/(6*np.pi))*alpha))/np.log(np.log2(n))-1/2 #IMPOTANT: MAYBE CHANGE TO LOG INSTEAD OF LOG2
        #and the error probability
        pr = 2*(1-norm.cdf(1/(4*t*alpha)*np.sqrt(6*np.pi/(r*(r+1)*m)))) #IMPORTANT: SHOULD BE UPDATED BASED ON NEW BOUND ON ROUNDING ERROR.
    if i==0:
        print("n    l    m      p       r   t  alpha    epsilon error Prob")
    print("%4d %4d %5d %7d %3d %3d %9.6f %5.3f %5.3f" %(n,l,m,p,r,t,alpha,epsilon,pr))
    
    
    #lwe = LWE(n,m,l,t,r,p,alpha)
    #keys = lwe.gen(elaborate)
    #priv,pub = keys
    #cipher = lwe.enc(v,pub,elaborate)
    #m = lwe.dec(cipher,priv,elaborate)