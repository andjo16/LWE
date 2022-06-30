'''
By Andreas H. Jørgensen (s202743)

This file is used to generate parameters for the LWE encryption scheme 
and generate empirical results on the decryption failure rate.
This file uses LWE64.py and FFPNG.py
'''

import sympy
from sympy.solvers import solve, nsolve
from sympy import Symbol, symbols, log, nonlinsolve, linsolve, Eq,N, Interval, Min,solveset,exp
from sympy.calculus.util import minimum
import math
from scipy.optimize import fmin, fsolve
import numpy as np
from LWE64 import LWE
import time
import os
import FFPNG
import multiprocessing as mp
import time

def findPandM(n,l,t,r):
    '''
    Function to find minimum p and m based on derived lower bounds
    '''
    def f(p):
        '''
        Function representing the intersection between the lower bound on p and m
        '''
        return (2*(n+l)*math.log(p[0],2)+200)/math.log(2*r+1,2)-p[0]/(2*t*r)
    
    p=fsolve(f,[100000]) #find p
    if(p>10):
            #p needs to be a prime
            p=math.ceil(p)
            #p=2000000 Manually setting P high
            if not sympy.isprime(p):
                p=sympy.nextprime(p)
            p=int(p)
            
            #Update m according to p
            m=mFunc(n,l,p,r) #math.ceil((2*(n+l)*math.log(p,2)+200)/math.log(2*r+1,2))
    else:
        print("WARNING: p is very small")
    return p,m

def mFunc(n,l,p,r):
    return math.ceil((2*(n+l)*math.log(p,2)+200)/math.log(2*r+1,2))

def findHermiteFactor(n,p,alpha):
    '''
    Used to find the minimum value for the hermite factor
    given alpha (and n and p)
    A hermite factor above 1.013 may give rise to an attack according to Gama and Nguyen
    '''
    return 2**((math.log(alpha/(1.5*math.sqrt(2*math.pi)),2)/(-2))**2/(n*math.log(p,2)))

def parameterInstation(n,t,r):
    '''
    Given a value of n, t and r this functions finds the
    remaining parameters according to the description in the report
    '''
    l=n
    p,m=findPandM(n,l,t,r) #first we find p and m

    #Now we find c (the exponent of gamma, as p>zeta*sqrt(ln(n)/n)
    c = math.log(p/(2*math.sqrt(math.log(n)/n)))/math.log(n)
    #c=math.floor(c * 100)/100.0
    
    #We can now find epsilon as the largest value that allows the upper and lower bounds on alpha to meet
    #This has been done by equalizing the two bounds and isolate epsilon
    epsilon = math.log(n**c*math.sqrt(math.log(n))/(n*t*math.sqrt((r*(r+1)*m)/(6*math.pi))))/(math.log(math.log(n)))-1/2
    #epsilon=math.floor(epsilon * 100)/100.0
    
    
    #Find the value of a that minimizes the upper bound for the probability of decryption failure
    def f(a):
        '''
        Decryption error as a function of a
        '''
        return l*2*math.exp(-2*m*r**2/a[0]**2)+l/((a[0]-1)/a[0]*(math.log(n)**(1/2+epsilon))/2)*math.exp(-(((a[0]-1)/a[0]*(math.log(n))**(1/2+epsilon))**2)/8)
    a = fmin(f,np.array([10]),disp=False)
    
    #Find remaining parameters
    alpha = 1/(t*math.sqrt((r*(r+1)*m/(6*math.pi)))*(math.log(n))**(1/2+epsilon))
    delta=findHermiteFactor(n,p,alpha)
    dim = math.ceil(math.sqrt(n*math.log(p)/math.log(delta)))
    keySize = math.ceil((n+l)*m*math.log(p))
    ciphertextBlowup = (1+n/l)*math.log(p)/math.log(t)

    #Return all results 
    ress = np.array([n,l,m,p,r,t,alpha,epsilon,c,a[0],delta,dim,keySize,ciphertextBlowup,f(a)*100])
    return ress


def parameterInstationV2(n,t,r,p):
    '''
    Given a value of n, t and r this functions finds the
    remaining parameters according to the description in the report
    '''
    l=n
    p2,m=findPandM(n,l,t,r) #first we find p and m
    if (p2>p):
        #Given p to small, we increase it to be p$
        p=p2
    else:
        if not sympy.isprime(p):
            p=sympy.nextprime(p)
        p=int(p)
        m=mFunc(n,l,p,r)
        
    #Now we find c (the exponent of gamma, as p>zeta*sqrt(ln(n)/n)
    c = math.log(p/(2*math.sqrt(math.log(n)/n)))/math.log(n)
    #c=math.floor(c * 100)/100.0
    
    #We can now find epsilon as the largest value that allows the upper and lower bounds on alpha to meet
    #This has been done by equalizing the two bounds and isolate epsilon
    epsilon = math.log(n**c*math.sqrt(math.log(n))/(n*t*math.sqrt((r*(r+1)*m)/(6*math.pi))))/(math.log(math.log(n)))-1/2
    #epsilon=math.floor(epsilon * 100)/100.0
    
    
    #Find the value of a that minimizes the upper bound for the probability of decryption failure
    def f(a):
        '''
        Decryption error as a function of a
        '''
        return l*2*math.exp(-2*m*r**2/a[0]**2)+l/((a[0]-1)/a[0]*(math.log(n)**(1/2+epsilon))/2)*math.exp(-(((a[0]-1)/a[0]*(math.log(n))**(1/2+epsilon))**2)/8)
    a = fmin(f,np.array([10]),disp=False)
    
    #Find remaining parameters
    alpha = 1/(t*math.sqrt((r*(r+1)*m/(6*math.pi)))*(math.log(n))**(1/2+epsilon))
    delta=findHermiteFactor(n,p,alpha)
    dim = math.ceil(math.sqrt(n*math.log(p)/math.log(delta)))
    keySize = math.ceil((n+l)*m*math.log(p))
    ciphertextBlowup = (1+n/l)*math.log(p)/math.log(t)

    #Return all results 
    ress = np.array([n,l,m,p,r,t,alpha,epsilon,c,a[0],delta,dim,keySize,ciphertextBlowup,f(a)*100])
    return ress

def onePercentageErrorParameters(parameters):
    '''
    This function finds the remaing parameters based on given n, r, t and epsilon
    The primary goal is to have parameter-sets with specific erorr probability.
    The epsilons's has been decided by visual inspection on a graph of the error probability.
    '''
    
    ress = np.empty((0,15))
    for params in parameters:
        n,r,t,epsilon = params
        l=n
    
        p,m = findPandM(n,l,t,r)
        
        #Now we find c (the exponent of gamma, as p>zeta*sqrt(ln(n)/n)
        c = math.log(p/(2*math.sqrt(math.log(n)/n)))/math.log(n)
        #c=math.floor(c * 100)/100.0
        
        
        #Find the value of a that minimizes the upper bound for the probability of decryption error
        def f(a):
            '''
            Decryption error as a function of a
            '''
            return l*2*math.exp(-2*m*r**2/a[0]**2)+l/((a[0]-1)/a[0]*(math.log(n)**(1/2+epsilon))/2)*math.exp(-(((a[0]-1)/a[0]*(math.log(n))**(1/2+epsilon))**2)/8)
        a = fmin(f,np.array([10]),disp=False)
        
        #Find remaining parameters
        alpha = 1/(t*math.sqrt((r*(r+1)*m/(6*math.pi)))*(math.log(n))**(1/2+epsilon))
        delta=findHermiteFactor(n,p,alpha)
        dim = math.ceil(math.sqrt(n*math.log(p)/math.log(delta)))
        keySize = math.ceil((n+l)*m*math.log(p))
        ciphertextBlowup = (1+n/l)*math.log(p)/math.log(t)
        
        ress = np.vstack([ress,[n,l,m,p,r,t,alpha,epsilon,c,a[0],delta,dim,keySize,ciphertextBlowup,f(a)*100]])
    return ress


def testDecryptErrorRateProcess(n,l,m,p,r,t,alpha,encCount,elaborate=False):
    '''
    This is used as a process for doing encCount encryptions and decryptions using the
    same key-pair which is generated by the function
    encCount: gives the number of messages that should be encrypted under the generated key
    '''
    lwe = LWE(n,m,l,t,r,p,alpha,elaborate)
    keys = lwe.gen(keyname=None,saveKey=False)
    priv,pub = keys
    letterError = 0
    msgError = 0
    for _ in range(encCount):
        #Encryption
        v = np.random.randint(0,t,l)
        if(elaborate):
            print("Message:", v)
        cipher = lwe.enc(v,pub)
        if cipher == None:
            print("ERROR")
            return
        if(elaborate):
            print("Ciphertext:",cipher)
        #Decryption
        mes = lwe.dec(cipher,priv)
        if m is None:
            print("ERROR")
            return
        mred = np.mod(mes,lwe.t)
        #Verification of results
        if(elaborate):
            print("Original message:")
            print(v)
            print("Decrypted message:")
            print(mes)
            print("Decrypted message modulo t:")
            print(mred)

        equal = v==mred
        errors = np.where(equal==False)
        errorsCount = len(errors[0])
        letterError+=errorsCount
        msgError += errorsCount>0
    return letterError,msgError

def testDecryptionErrorRateMultiprocessing(n,l,m,p,r,t,alpha,repetitions,encCount,generateKey=False,elaborate=False):
    '''
    Function for starting a thread pool, and start a new thread for each 'repetitions'
    Each of the repetitions threads uses the function 'testDecryptErrorRateProcess()' to do encCount encryptions for
    freshly generated key-pair
    repetitions: gives the number of times a new key should be generated
    encCount: gives the number of messages that should be encrypted under each key
    '''
    
    timestr = time.strftime("%y%m%d-%H,%M,%S")
    tmpFileName = "tmpDecError" + timestr + ".txt"
    resultFileName = "DecErrorProb"+str(n)+","+str(m)+".txt"
    
    if not os.path.exists('results/'):
        os.makedirs('results')
    tmpFilePath = os.path.join('results', tmpFileName)
    resultFilePath = os.path.join('results', resultFileName)
    
    with mp.Pool(processes =mp.cpu_count()+1) as pool:
        processes = list()
        results = list()
        for _ in range(repetitions):
            res = pool.apply_async(testDecryptErrorRateProcess,(n,l,m,p,r,t,alpha,encCount,False))
            results.append(res)

        letterError = 0 #total number of letters having an error
        letterCount = 0
        msgError = 0 #Total number of messages in which at least one letter is errored
        msgCount = 0
        for i,res in enumerate(results):
            lErr,mErr = res.get()
            letterError += lErr
            letterCount = (l*(i+1))*encCount
            letterErrorFraction = letterError/letterCount
            msgError += mErr
            msgCount = (i+1)*encCount
            msgErrorFraction = msgError/msgCount
            if (i+1)%10==0:
                print("Letter error: %d of %d: %.2e%% || Message error: %d of %d: %.2e%%" %(letterError, letterCount, letterErrorFraction*100, msgError,msgCount,msgErrorFraction*100))
                with open(tmpFilePath,"a") as f:
                    f.write("Letter Error: %4d of %6d: %.2e%% || Message error: %d of %d: %.2e%% - (n,l,m,p,r,t,alpha): (%d, %d, %d, %d, %d, %d, %8.4f)\n" %(letterError, letterCount, letterErrorFraction*100,msgError,msgCount,msgErrorFraction*100,n,l,m,p,r,t,alpha))
        with open(resultFilePath,"a") as f:
            f.write("Letter error: %d of %d%% || Message error: %d of %d: %.2e%%: %.2e\n %d, %d, %d, %d, %d, %d, %8.4f\n" %(letterError, letterCount, letterErrorFraction*100,msgError,msgCount,msgErrorFraction*100,n,l,m,p,r,t,alpha))

    return letterErrorFraction*100,msgErrorFraction*100


def empiricalError(results,repetitions,encCount,generateKey,elaborate):
    '''
    This function is used to get an empirical estimate for the probability
    of deccryption failure for each set of parameters given in 'results'
    Estimates is done by generating 'repetitions' many keys and then
    do 'encCount' encryptions pr. key.
    Both the fraction of the total number of letter encrypted and the total number
    of messages (consisting of n letters) giving decryption failures are found
    repetitions: the number of keys to be generated
    encCount: the number of messages to be generated under each key
    '''
    experimentalError = [] #Expected probability of encryption error
    forcedError = [] #Probability of being able to enforce encryption error using a LWE-oracle
    for parameters in results:
        n,l,m,p,r,t,alpha = parameters[0:7]
        n=int(n);l=int(l);m=int(m);p=int(p);r=int(r);t=int(t)        
        letterFraction,msgFraction = testDecryptionErrorRateMultiprocessing(n,l,m,p,r,t,alpha,repetitions,encCount,generateKey,elaborate)
        experimentalError.append([letterFraction,msgFraction,repetitions*encCount,encCount])
    return experimentalError

def forcedError(results,repetitions,generateKey,elaborate):
    '''
    The empirical probability of being able to enforce a decryption failure
    given an LWE-oracle (that is, given acccess to the error vector E).
    The estimates is done based on 'repetitions' number of letters for each set
    of parameters in 'results'
    repetitions: The number of letters used to estitmate the success probability
    '''
    experimentalError = [] #Probability of being able to enforce encryption error using a LWE-oracle
    for parameters in results:
        n,l,m,p,r,t,alpha = parameters[0:7]
        n=int(n);l=int(l);m=int(m);p=int(p);r=int(r);t=int(t)
        experimentalError.append([FFPNG.errorProbParallel(m,p,t,r,alpha,repetitions,elaborate),repetitions])
    return experimentalError

def mergeColumns(combinedTitles,results1,titles1, results2=None, titles2=None):
    '''
    This function should be used to merge the results of two numpy columns/vectors
    The data is merged based on whether the corresponding entries in 'titles1' and 'titles2' are the same
    Data that is available in both vectors, but doesn't match will be updated with the data in results1. A warning will be printed in that case
    combinedTitles: the combined set of titles. This determines the order of the data in the combined column
    If results2 and titles2 = None, then we just need to extend results1 to match Combined Titles. Unknown values will be given value NaN
    '''
    combined = np.empty((len(combinedTitles))) #merged vector
    if results1 is None or titles2 is None:
        #we just need to extend the columns to match combinedTitles
        for i,t in enumerate(combinedTitles):
            if t in titles1:
                #add the value
                j = np.where(titles1 == t)[0][0]
                combined[i]=results1[j]
            else:
                #add NaN
                combined[i]=float("nan")
    else:
        #We need to combine two columns
        for i,t in enumerate(combinedTitles):
            if t in titles1 and t in titles2:
                #entry in both vectors: merge needed
                j1 = np.where(titles1 == t)[0][0]
                j2 = np.where(titles2 == t)[0][0]
                if(not math.isclose(results1[j1],results2[j2])):
                    #If not equal, print warning and use value in results1
                    print("Warning: '%s' of two vectors is not equal: %s!=%s. The first one is used" %(t,results1[j1],results2[j2]))
                combined[i]=results1[j1]
            elif t in titles1:
                #entry only in first vector: add to result
                j = np.where(titles1 == t)[0][0]
                combined[i]=results1[j]
            elif t in titles2:
                #entry only in second vector: add to result
                j = np.where(titles2 == t)[0][0]
                combined[i]=results2[j]
    return combined

def mergeResultsWithFile(results,titles,inputFile):
    '''
    This function will take a numpy-matrix of results, and a numpy-matrix from inputFile
    of previous results, and combine the two matrices.
    The resulting matrix is returned
    A column from the results-matrix is merged with a column from the inputFile matrix in case n, m, p and alpha matches
    '''
    print("###Started merging###")
    loadedResults,loadedTitles = loadResults(inputFile)
    
        
    #merge the two titles arrays. The order in loadTitles takes precedence
    combinedTitles = loadedTitles
    tmp = [t for t in titles if t not in combinedTitles]
    combinedTitles = np.append(combinedTitles,tmp)
    
    #Find the indices in result- and loadedResult-matrix where the values for n, m, p and alpha are found
    cmpParameters = ["$n$", "$m$", "$p$", "$\\alpha$"]
    cmpIndices = []
    for param in cmpParameters:
        cmpIndices.append((np.where(titles == param)[0][0],np.where(loadedTitles == param)[0][0]))

    #Find matching columns
    matches = [] #List of tuples indicating which columns of results that matches which columns of loadedResults
    resSub = results[list(zip(*cmpIndices))[0],:] #extract rows of interest for comparison
    loadedSub = loadedResults[list(zip(*cmpIndices))[1],:] #extract rows of interest for comparison

    for i,resCol in enumerate(resSub.T):
        #For each column in results, determine if it matches any of the columns in loadedResults
        resCol=resCol.reshape(-1,1)
        tmpList = np.where(np.all(loadedSub == resCol,axis=0)==True)[0]
        if len(tmpList)>0:
            if len(tmpList)>1:
                print("ERROR too many matched columns")
                print("Aborts: Hopefully backup was implemented")
                return
            matches.append((i,tmpList[0]))
    
    combinedResults = np.zeros((len(combinedTitles),0))
    if len(matches)>0:
        resMatchInd = np.array(list(zip(*matches))[0])
        loadMatchInd = np.array(list(zip(*matches))[1])
    else:
        resMatchInd = []
        loadMatchInd = []
    #We want to keep the order of the columns in the loadedResults-matrix and then append
    #new columns to the end of the matrix
    for i,loadCol in enumerate(loadedResults.T):
        #For each column in the loadedResult-matrix
        if i in loadMatchInd:
            #the same column is in both matrices
            #We will combine the columns and add the resulting column to the resulting matrix
            j = resMatchInd[np.where(loadMatchInd == i)[0]][0]
            print("Column %d of old data matches column %d of new data" % (i,j))
            col = mergeColumns(combinedTitles,results.T[j],titles,loadedResults.T[i],loadedTitles) #the new data is prefered in case of mismatch in data
            combinedResults = np.hstack((combinedResults,col.reshape(-1,1)))
        else:
            #if a column is only in the loaded matrix
            print("Column %d of old data is unique" %i)
            col = mergeColumns(combinedTitles,loadedResults.T[i],loadedTitles)
            combinedResults = np.hstack((combinedResults,col.reshape(-1,1)))
    
    #now we add all columns of new data not already merged
    ind = [i for i in range(results.shape[1]) if i not in resMatchInd]
    for i in ind:
        print("Column %d of new data is unique" %i)
        col = mergeColumns(combinedTitles,results.T[i],titles)
        combinedResults = np.hstack((combinedResults,col.reshape(-1,1)))
    
    return (combinedResults,combinedTitles)

def saveResults(res,titles,filePath):
    '''
    Saves the titles array and result matrix in the file filePath
    '''
    with open(filePath,"w+b") as f:
        np.save(f,titles)
        np.save(f,res)

def loadResults(filePath):
    '''
    Loads the titles array and result matrix from the file filePath
    '''
    with open(filePath,"r+b") as f:
        titles = np.load(f)
        res = np.load(f)
    return (res,titles)
    
def matrixToLatex(m,titles,formats):
    '''
    Converts a matrix 'm' of results to a latex table for report
    '''
    string = ""
    rows,columns = m.shape
    string += "\\begin{tabular}{|"
    for _ in range(columns+2):
        string+= "c|"
    string+= "}\n"
    string+="\hline\n"
    string+="\multirow{8}{*}{Parameters}\n"
    for i,row in enumerate(m):
        if i<=9:
            tmp = "& %15s " %titles[i]#%row[0]
        else:
            tmp = "\\multicolumn{2}{|c|}{\\makecell{%15s}}" %titles[i]
        for d in row:
            if(np.isnan(d)):
                tmp += " & %f" %d
            else:
                if "prob" in titles[i]:
                    if (d<0.001):
                        tmp += " & %8.2e" %d
                    else:
                        tmp += " & %8.4f" %d
                    tmp += "\\%"
                else:
                    tmp += " & "+formats[titles[i]] %d            
        if i ==6 or i>=9:
            tmp += "\\\\\n\hline\n"
            if i==6:
                tmp +="\multirow{3}{*}{\\shortstack{Meta\\\\ parameters}}\n"
        else:
            tmp+= "\\\\\n\cline{2-%d}\n" %(columns+2)
        if i<len(titles)-2 and ("FFP-NG" in titles[i+1] or "letter" in titles[i+1]):
            tmp+="\hline\n"
        string += tmp
    string+="\\end{tabular}"
    return string


    


def main():
    '''
    Encryption times (empirical decryption error estimation)
    n=512: 10000 reps * 1 key reuse: 150 hours
           1000 reps * 10 key reuse: 15 hours 
           100 reps * 100 key reuse: 1.67 hours
    n=384: 10000 reps * 1 key reuse: 50 hours
           1000 reps * 10 key reuse: 3,33 hours
           100 reps * 100 key reuse: 0.5 hours
    n=256: 10000 reps * 1 key reuse: 4.15 hours
    n=128: 10000 reps * 1 key reuse: 0.5 hours
    
    FFP-NG times (empirical forced decryption error estimation)
    n=128: 1,000,000 reps: 20 minutes
    n=256: 1,000,000 reps: 25 minutes
    n=384: 1,000,000 reps: 40-50 minutes
    n=512: 1,000,000 reps: 50 minutes
    '''
    
    
    
    t=2
    r=1
    '''
    n=64
    res = parameterInstation(n,t,r)
    results =res #np.vstack([results,res])
    
    n=100
    res = parameterInstation(n,t,r)
    results =np.vstack([results,res])
    '''
    n=128
    res = parameterInstation(n,t,r)
    #results =res#np.vstack([results,res])
    
    n=256
    res=parameterInstation(n,t,r)
    #results = np.vstack([results,res])
    
    n=384
    res=parameterInstation(n,t,r)
    #results = np.vstack([results,res])
    
    n=512
    res=parameterInstation(n,t,r)
    #results = np.vstack([results,res])
   
    n=256
    res=parameterInstation(n,2,50)
    results =res
    res=parameterInstation(n,50,1)
    results = np.vstack([results,res])
    res = parameterInstationV2(n,t,r,200000) #200k
    results = np.vstack([results,res])
    res = parameterInstationV2(n,2,50,2000000) #2mio
    results = np.vstack([results,res])
    res = parameterInstationV2(n,50,1,2000000) #2mio
    results = np.vstack([results,res])
    
    #Parameters to achieve a specicif decryption error rate
                             #n,r,t,epsilon
    parameters = np.array([ [256,1,2,0.82], #0.1%
                            [256,1,2,0.75], #1%
                            [256,1,2,0.67], #10%
                            [128,1,2,0.92], #0.1%
                            [128,1,2,0.84], #1%
                            [128,1,2,0.74], #10%
                            [384,1,2,0.77], #0.1%
                            [384,1,2,0.71], #1%
                            [384,1,2,0.64], #10%
                            [512,1,2,0.75], #0.1%
                            [512,1,2,0.69], #1%
                            [512,1,2,0.62]  #10%
                            ])
    #Epsilon is chosen by visual inspection of graphs

    #results =np.vstack([results,onePercentageErrorParameters(parameters)])

    titles = ["$n$","$\\ell$","$m$","$p$","$r$","$t$","$\\alpha$","$\\epsilon$","$c$","$a$","$\\delta$","Lattice dim. \\\\ in attack","Public key size","Blowup factor","Error prob. \\\\ upper bound"]
    formats = ["%8d","%8d", "%8d","%8d","%8d","%8d","%8.4f","%8.2f","%8.2f","%8.2f","%8.4f","%5d","%8.2e","%8.2f","%8.3e"]
    titles = np.array(titles)
    titlesDict = dict(zip(titles,formats))
    titlesDict.update(dict(zip(["Empirical letter \\\\ error prob.","Empirical message \\\\ error prob.","Encryption Count","Key reuse"],["%8.3e","%8.3e","%d","%d"])))
    titlesDict.update(dict(zip(["Empirical FFP-NG \\\\ error prob.","Sample Count"],["%8.3e","%d"])))
    titlesDict.update(dict(zip(["Prob. ratio"],["%8.2f"])))
    titlesDict.update(dict(zip(["gen", "enc", "dec","size","blowup"],["%8.2e","%8.2e","%8.2e","%8.2e","%8.2f"])))
    
    elaborate=False
    generateKey=True
    
    findErrorProb = False #estimate decryption error probability
    findForcedErrorProb = False #estimate success probability of forcing decryption given LWE-oracle

    if(results.ndim == 1):
        results = np.reshape(results,(-1,results.size))

    if findErrorProb:
        repetitions=10000 #number of keys #takes about 90 minutes with 10.000 without repeated key generation
        encryptionsPrKey = 1 #number of messages to be encrypted for each key
        expError = empiricalError(results,repetitions,encryptionsPrKey,generateKey,elaborate)
        results=np.hstack((results,np.array(expError).reshape(-1,4)))
        titles=np.append(titles, ["Empirical letter \\\\ error prob.","Empirical message \\\\ error prob.","Encryption Count","Key reuse"])
        formats.append("%8.3e")
        formats.append("%8.3e")
        formats.append("%d")
        formats.append("%d")
    
    if findForcedErrorProb:
        repetitions=1000000
        err = forcedError(results,repetitions,generateKey,elaborate)
        results=np.hstack((results,np.array(err).reshape(-1,2)))
        titles = np.append(titles,["Empirical FFP-NG \\\\ error prob.","Sample Count"])
        formats.append("%8.3e")
        formats.append("%d")
    
    results=results.T
    
    #Prepare storing results
    timestr = time.strftime("%y%m%d-%H,%M,%S")
    tmpResultsFile = "tmpResults" + timestr + ".npy" #This is a temporary file used for results from this run only
    resultFile = "combinedResults.npy" #File where results from every run is combined
    resultFileFinal = "combinedResultsFinal.npy" #Final file where results from every run is found - this file is not written to
    tmpLatexFile = "tmpLatex - " + timestr + ".txt" #Latex table version of result from this run only
    latexFile = "combinedLatex.txt" #Latex table version of results from every run combined
    
    folder = "results2/"
    if not os.path.exists(folder):
        os.makedirs(folder)
    tmpResultsPath = os.path.join(folder,tmpResultsFile)
    resultPath = os.path.join(folder,resultFile) #combined results
    resultPathFinal = resultFileFinal #File containing all presented data. Only to be read from.
    tmpLatexPath = os.path.join(folder,tmpLatexFile)
    latexPath = os.path.join(folder,latexFile) #combiend results
    
    #Store temporary results before combinining
    saveResults(results,titles,tmpResultsPath)
    table = matrixToLatex(results,titles,titlesDict)
    print(table)
    with open(tmpLatexPath,"w") as f:
        f.write(table)
    
    #Combine results
    if(os.path.exists(resultPathFinal)):
        results,titles = mergeResultsWithFile(results,titles,resultPathFinal)
    
    #Used to find the ratio between theoretical upper bound and empirical decryption failure rate
    if("Error prob. \\\\ upper bound" in titles and "Empirical message \\\\ error prob." in titles):
        i = np.where(titles=="Error prob. \\\\ upper bound")[0][0]
        j = np.where(titles=="Empirical message \\\\ error prob.")[0][0]
        tmp = results[[i,j],:]
        probRatio = results[i,:]/results[j,:]
        if("Prob. ratio" in titles):
            #if ratio already in table, then recalculate and insert at found location
            i = np.where(titles=="Prob. ratio")[0][0]
            results[i,:]=probRatio #recalculate ratio and insert
        else:
            #if ratio not already in table, then insert it at correct location
            print("#####################Prob ratio not found")
            titles = np.append(titles,["Prob. ratio"])
            formats.append("%8.2f")
            results=np.vstack((results,np.array(probRatio).reshape(1,-1)))
            results[[19,20,21],:]=results[[21,19,20],:]
            titles[[19,20,21]]=titles[[21,19,20]]
    
    #results[:,[4,5,6,7,8,9]]=results[:,[7,8,9,4,5,6]]
    results=results[:,[1,7,8,9,16,17,18,19,20]]
    
    n = results[0]
    print(n)
    l = results[1]
    print(l)
    m = results[2]
    p = results[3]
    r = results[4]
    t = results[5]
    
    gen = n*m*l*(np.log(p)**2)
    enc = (n+l)*m*(np.log(p)**2)
    dec = n*l*(np.log(p)**2)
    #size = (n+l)*m*(np.log(p))
    #blowup = (1+n/l)*(np.log(p)/np.log(t))
    print(type(titles))
    titles = np.concatenate((titles,np.array(["gen", "enc", "dec"])),axis=None)
    #titles = np.concatenate((titles,np.array(["gen", "enc", "dec","size","blowup"])),axis=None)
    #formats= np.concatenate((formats,np.array(["%8.2f","%8.2f","%8.2f","%8.2f","%8.2f"])), axis=None)
    results=np.vstack((results,np.array(gen).reshape(1,-1)))
    results=np.vstack((results,np.array(enc).reshape(1,-1)))
    results=np.vstack((results,np.array(dec).reshape(1,-1)))
    #results=np.vstack((results,np.array(size).reshape(1,-1)))
    #results=np.vstack((results,np.array(blowup).reshape(1,-1)))
    
    results = results[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,22,23,24],:]
    titles = titles[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,22,23,24]]
    
    #Store combined results
    saveResults(results,titles,resultPath)
    #results=results[:,[0,1,2,3,7,8,9]] #For table in report
    #results=results[:,4:] #For table in appendix
    table = matrixToLatex(results,titles,titlesDict)
    print(table)
    with open(latexPath,"w") as f:
        f.write(table)
    
if __name__ == "__main__":
    main()
