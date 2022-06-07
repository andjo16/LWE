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
#x = Symbol('x')
#print(solve(x**2 - 1, x))


'''
This file uses LWE64.py and FFPNG.py
'''


def findPandM(n,l,t,r):
    '''
    Function to find minimum p and m based on derived lower bounds
    '''
    def f(p):
        '''
        Function representing the intersection between the lower bound on p and m
        '''
        return (2*(n+l)*math.log(p[0],2)+200)/math.log(2*r+1,2)-p[0]/(2*t*r)
    
    p=fsolve(f,[10000]) #find p
    if(p>10):
            #p needs to be a prime
            p=math.ceil(p)
            if not sympy.isprime(p):
                p=sympy.nextprime(p)
            p=int(p)
            #Update m according to p
            #m=math.ceil(3*(n+l)*math.log(p)/math.log(2*r+1)) #Remember to update m
            m=math.ceil((2*(n+l)*math.log(p,2)+200)/math.log(2*r+1,2))
    else:
        print("WARNING: p is very small")
    return p,m

def findHermiteFactor(n,p,alpha):
    '''
    Used to find the minimum value for the hermite factor
    given alpha (and n and p)
    A hermite factor above 1.013 may give rise to an attack according to Gama and Nguyen
    '''
    return 2**((math.log(alpha/(1.5*math.sqrt(2*math.pi)),2)/(-2))**2/(n*math.log(p,2)))

def parameterInstation(n=128,t=2,r=1):
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
    #print(n,p,alpha,delta,dim)
    keySize = math.ceil((n+l)*m*math.log(p))
    ciphertextBlowup = (1+n/l)*math.log(p)/math.log(t)

    #Return all results 
    ress = np.array([n,l,m,p,r,t,alpha,epsilon,c,a[0],delta,dim,keySize,ciphertextBlowup,f(a)*100])
    return ress

def onePercentageErrorParameters(parameters):
    '''
    This function finds the remaing parameters based on given n, r, t and epsilon
    The primary goal is to have parameter-sets with specific erorr probability.
    The epsilons's has been decided by visual inspection of the error probability.
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
        #print(n,p,alpha,delta,dim)
        keySize = math.ceil((n+l)*m*math.log(p))
        ciphertextBlowup = (1+n/l)*math.log(p)/math.log(t)
        
        ress = np.vstack([ress,[n,l,m,p,r,t,alpha,epsilon,c,a[0],delta,dim,keySize,ciphertextBlowup,f(a)*100]])
    return ress


def testDecryptErrorRateProcess(n,l,m,p,r,t,alpha,encCount,elaborate=False):
    '''
    encCount gives the number of messages that should be encrypted under the generated key
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
    repetitions gives the number of times a new key should be generated
    encCount gives the number of messages that should be encrypted under each key
    '''
    
    timestr = time.strftime("%y%m%d-%H,%M,%S")
    tmpFileName = "tmpDecError" + timestr + ".txt"
    resultFileName = "DecErrorProb"+str(n)+","+str(m)+".txt"
    
    if not os.path.exists('results/'):
        os.makedirs('results')
    tmpFilePath = os.path.join('results', tmpFileName)
    resultFilePath = os.path.join('results', resultFileName)
    
    #print(mp.cpu_count())
    with mp.Pool(processes =mp.cpu_count()+1) as pool:
        #dataQueue = None #mp.Queue(100)
        processes = list()
        results = list()
        for _ in range(repetitions):
            res = pool.apply_async(testDecryptErrorRateProcess,(n,l,m,p,r,t,alpha,encCount,False))
            results.append(res)
            #proc = mp.Process(target=testDecryptErrorRateProcess, args=(dataQueue,n,l,m,p,r,t,alpha,False))
            #processes.append(proc)
            #proc.start()

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

    #errorFraction = sum/(l*repetitions)
    return letterErrorFraction*100,msgErrorFraction*100
    
    
    
    #for proc in processes:
    #    proc.join()
    
    #for _ in range(repetitions):
    #    res = dataQueue.get()
    #    print(res)
    
    
    
    # #Seems to work, but creates a new process for each task
    # dataQueue = mp.Queue(100)
    # processes = list()
    # for _ in range(repetitions):
        # proc = mp.Process(target=testDecryptErrorRateProcess, args=(dataQueue,n,l,m,p,r,t,alpha,False))
        # processes.append(proc)
        # proc.start()
    
    # #for proc in processes:
    # #    proc.join()
    
    # for _ in range(repetitions):
        # res = dataQueue.get()
        # print(res)
    
    
    
    
    
    # #One solution that works, but only gives all results in the end
    # pool = mp.Pool(processes =4)
    # results = pool.starmap(testDecryptErrorRateProcess,[(n,l,m,p,r,t,alpha,False)]*repetitions)
    # print(results)
    
    """
    totalError = 0
    
    for i in range(repetitions):
        
        totalError+=len(errors[0])
        totalErrorFraction = totalError/(l*(i+1))
        if (i+1)%1==0:
            print("Error: %d of %d: %.2e%%" %(totalError, l*(i+1), totalErrorFraction*100))
            with open(tmpFilePath,"a") as f:
                f.write("Error: %4d of %6d: %.2e%% - (n,l,m,p,r,t,alpha): (%d, %d, %d, %d, %d, %d, %8.4f)\n" %(totalError,l*(i+1), totalErrorFraction*100,n,l,m,p,r,t,alpha))
    with open(resultFilePath,"a") as f:
        f.write("Error: %d of %d%%: %.2e\n %d, %d, %d, %d, %d, %d, %8.4f\n" %(totalError,l*(i+1), totalErrorFraction*100,n,l,m,p,r,t,alpha))
    return totalErrorFraction*100
    """

def testDecryptionErrorRate(n,l,m,p,r,t,alpha,repetitions,generateKey=False,elaborate=False):
    timestr = time.strftime("%y%m%d-%H,%M,%S")
    tmpFileName = "tmpDecError" + timestr + ".txt"
    resultFileName = "DecErrorProb"+str(n)+","+str(m)+".txt"
    
    if not os.path.exists('results/'):
        os.makedirs('results')
    tmpFilePath = os.path.join('results', tmpFileName)
    resultFilePath = os.path.join('results', resultFileName)

    ##Preparing/setup
    keyname="key"+str(n)+","+str(m)+".npy" #name of file to save/load the generated key
    lwe = LWE(n,m,l,t,r,p,alpha,elaborate)
    
    #DO I NEED TO GENERATE A NEW KEY EACH TIME????
    if not os.path.exists(keyname):
    #if generateKey:
        keys = lwe.gen(keyname)
        priv,pub = keys
    else:
        with open(keyname,"rb") as f:
            priv = np.load(f)
            pub = np.load(f,allow_pickle=True),np.load(f,allow_pickle=True)
    
    totalError = 0
    for i in range(repetitions):
        #DO I NEED TO GENERATE A NEW KEY EACH TIME????
        
        #Generate key
        keys = lwe.gen(keyname,saveKey=False)
        priv,pub = keys
        
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
        if(elaborate):
            print("Found %d errors" %len(errors[0]))
            if(len(errors[0])<20):
                for i in errors[0]:
                    print("Error at index %d: %d != %d" %(i,v[i],mred[i]))
        totalError+=len(errors[0])
        totalErrorFraction = totalError/(l*(i+1))
        if (i+1)%1==0:
            print("Error: %d of %d: %.2e%%" %(totalError, l*(i+1), totalErrorFraction*100))
            with open(tmpFilePath,"a") as f:
                f.write("Error: %4d of %6d: %.2e%% - (n,l,m,p,r,t,alpha): (%d, %d, %d, %d, %d, %d, %8.4f)\n" %(totalError,l*(i+1), totalErrorFraction*100,n,l,m,p,r,t,alpha))
    with open(resultFilePath,"a") as f:
        f.write("Error: %d of %d%%: %.2e\n %d, %d, %d, %d, %d, %d, %8.4f\n" %(totalError,l*(i+1), totalErrorFraction*100,n,l,m,p,r,t,alpha))
    return totalErrorFraction*100
    




def empiricalError(results,repetitions,encCount,generateKey,elaborate):
    '''
    The empirical probability of getting af row with an error.
    repetitions is the number of keys to be generated
    encCount is the number of messages to be generated under each key
    '''
    experimentalError = [] #Expected probability of encryption error
    forcedError = [] #Probability of being able to enforce encryption error using a LWE-oracle
    for parameters in results:
    #res = results[0]
    #for res in results[1:]:
        n,l,m,p,r,t,alpha = parameters[0:7]
        n=int(n);l=int(l);m=int(m);p=int(p);r=int(r);t=int(t)
        #repetitions=10 #takes about 90 minutes with 10.000 without repeated key generation
        
        #errorFraction = testDecryptionErrorRate(n,l,m,p,r,t,alpha,repetitions,generateKey,elaborate)
        letterFraction,msgFraction = testDecryptionErrorRateMultiprocessing(n,l,m,p,r,t,alpha,repetitions,encCount,generateKey,elaborate)
        experimentalError.append([letterFraction,msgFraction,repetitions*encCount,encCount])
        #print(errorFraction)
    return experimentalError

def forcedError(results,repetitions,generateKey,elaborate):
    '''
    The empirical probability of being able to enforce a decryption error
    given an LWE-oracle
    '''
    experimentalError = [] #Probability of being able to enforce encryption error using a LWE-oracle
    for parameters in results:
        n,l,m,p,r,t,alpha = parameters[0:7]
        n=int(n);l=int(l);m=int(m);p=int(p);r=int(r);t=int(t)
        #repetitions=10000
        #t = time.time()
        experimentalError.append([FFPNG.errorProbParallel(m,p,t,r,alpha,repetitions,elaborate),repetitions])
        #print("time: %0.2f" %(time.time()-t))
    return experimentalError

def mergeColumns(combinedTitles,results1,titles1, results2=None, titles2=None):
    '''
    This function should be used to merge the results of two numpy columns/vectors
    Data that is available in both vectors, but doesn't match will be updated with the data in results 1. A warning will be printed in that case
    combinedTitles is the combined set of titles
    If results2 and titles2 = None, then we just need to extend results1 to match Combined Titles. Unknown values will be given NaN (float("nan"))
    '''
    combined = np.empty((len(combinedTitles)))
    if results1 is None or titles2 is None:
        #we just need to extend the columns to match combinedTitles
        for i,t in enumerate(combinedTitles):
            if t in titles1:
                #add the value
                j = np.where(titles1 == t)[0][0]
                combined[i]=results1[j]
            else:
                #add NaN
                #j = np.where(titles1 == t)[0][0]
                combined[i]=float("nan")
    else:
        #We need to combine two columns
        for i,t in enumerate(combinedTitles):
            if t in titles1 and t in titles2:
                j1 = np.where(titles1 == t)[0][0]
                j2 = np.where(titles2 == t)[0][0]
                if(not math.isclose(results1[j1],results2[j2])):
                    #If not equal, print warning and use results1
                    print("Warning: '%s' of two vectors is not equal: %s!=%s. The first one is used" %(t,results1[j1],results2[j2]))
                combined[i]=results1[j1]
            elif t in titles1:
                #add the value
                j = np.where(titles1 == t)[0][0]
                combined[i]=results1[j]
            elif t in titles2:
                #add the value
                j = np.where(titles2 == t)[0][0]
                combined[i]=results2[j]
    return combined

def mergeResultsWithFile(results,titles,inputFile):
    '''
    This function will take a np-matrix of results, and a np-matrix from inputFile
    of previous results, and combine the two matrices.
    The resulting matrix is returned
    A column from the results-matri is merged with a column from the inputFile matrix in case n, m, p and alpha matches
    '''
    print("###Started merging###")
    loadedResults,loadedTitles = loadResults(inputFile)
    
        
    #merge the two titles arrays. The order in loadTitles takes precedence
    combinedTitles = loadedTitles
    tmp = [t for t in titles if t not in combinedTitles]
    combinedTitles = np.append(combinedTitles,tmp)
    
    #Find the indices in result and loadedResult-matrix where the values for n, m, p and alpha are found
    cmpParameters = ["$n$", "$m$", "$p$", "$\\alpha$"]
    cmpIndices = []
    for param in cmpParameters:
        cmpIndices.append((np.where(titles == param)[0][0],np.where(loadedTitles == param)[0][0]))

    #Find matching columns
    matches = [] #List of tuples indicating which columns of results that matches which columns of loadedResults
    resSub = results[list(zip(*cmpIndices))[0],:] #extract rows of interest for comparison
    loadedSub = loadedResults[list(zip(*cmpIndices))[1],:] #extract rows of interest for comparison

    for i,resCol in enumerate(resSub.T):
        resCol=resCol.reshape(-1,1)

        tmpList = np.where(np.all(loadedSub == resCol,axis=0)==True)[0]#np.where(loadedSub == resCol)[0]
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
    #We want to keep the order of the columns in the loaded data-matrix and then append
    #new columns to the end of the matrix
    for i,loadCol in enumerate(loadedResults.T):
        if i in loadMatchInd:
            #If the same column is in both matrices
            j = resMatchInd[np.where(loadMatchInd == i)[0]][0]
            print("Column %d of old data matches column %d of new data" % (i,j))
            col = mergeColumns(combinedTitles,results.T[j],titles,loadedResults.T[i],loadedTitles) #the new data is prefered in case of mismatch in data
            combinedResults = np.hstack((combinedResults,col.reshape(-1,1)))
        else:
            #if a column is only in the loaded matrix
            print("Column %d of old data is unique" %i)
            col = mergeColumns(combinedTitles,loadedResults.T[i],loadedTitles)
            combinedResults = np.hstack((combinedResults,col.reshape(-1,1)))
    
    ind = [i for i in range(results.shape[1]) if i not in resMatchInd] #all columns of new data not already merged
    for i in ind:
        #append every column in the new result matrix
        print("Column %d of new data is unique" %i)
        col = mergeColumns(combinedTitles,results.T[i],titles)
        combinedResults = np.hstack((combinedResults,col.reshape(-1,1)))
    
    return (combinedResults,combinedTitles)

def saveResults(res,titles,filePath):
    '''
    Saves the titles array and result matrix in the file fileName
    '''
    #if not os.path.exists('results/'):
    #    os.makedirs('results')
    #filePath = os.path.join('results', fileName)
    with open(filePath,"w+b") as f:
        np.save(f,titles)
        np.save(f,res)

def loadResults(filePath):
    #if not os.path.exists('results/'):
    #    os.makedirs('results')
    #filePath = os.path.join('results', fileName)
    with open(filePath,"r+b") as f:
        titles = np.load(f)
        res = np.load(f)
    return (res,titles)
    
def matrixToLatex(m,titles,formats):
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
            #if i in [12,13,14,15,16]:
            
        if i ==6 or i>=9:
            tmp += "\\\\\n\hline\n"
            if i==6:
                tmp +="\multirow{3}{*}{\\shortstack{Meta\\\\ parameters}}\n"
            #if i==9:
            #    tmp +="\multirow{3}{*}{}\n"
        else:
            tmp+= "\\\\\n\cline{2-%d}\n" %(columns+2)
        #if i == (len(titles)-2):
        if i<len(titles)-2 and "FFP-NG" in titles[i+1]:
            tmp+="\hline\n"
        string += tmp
        
    string+="\\end{tabular}"
    #print(string)
    return string


def matrixToLatexOLD(m,titles,formats):
    string = ""
    rows,columns = m.shape
    string += "\\begin{tabular}{|"
    for _ in range(columns+2):
        string+= "c|"
    string+= "}\n"
    string+="\hline\n"
    string+="\multirow{8}{*}{Parameters}\n"
    for i,row in enumerate(m):
        print(i)
        tmp = "& %15s " %titles[i]#%row[0]
        for d in row:
            tmp += " & "+formats[i] %d
            if i in [12,13,14,15]:
                tmp += "\\%"
        if i in [6,9,rows-1]:
            tmp += "\\\\\n\hline\n"
            if i==6:
                tmp +="\multirow{3}{*}{\\shortstack{Meta\\\\ parameters}}\n"
            if i==9:
                tmp +="\multirow{3}{*}{}\n"
        else:
            tmp+= "\\\\\n\cline{2-%d}\n" %(columns+2)
        string += tmp
    string+="\\end{tabular}"
    #print(string)
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
    
    
    
    #results = np.array(["$n$","$\ell$","$m$","$p$","$r$","$t$","$\\alpha$","$\epsilon$","$c$","$a$","Public key size","Blowup factor","Error prob."])
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
    
    #Parameters to achieve a specicif decryption error rate
                             #n,r,t,epsilon
    parameters = np.array([ #[256,1,2,0.82], #0.1%
                            #[256,1,2,0.75], #1%
                            [256,1,2,0.67], #10%
                            #[128,1,2,0.92], #0.1%
                            #[128,1,2,0.84], #1%
                            [128,1,2,0.74], #10%
                            #[384,1,2,0.77], #0.1%
                            #[384,1,2,0.71], #1%
                            [384,1,2,0.64], #10%
                            #[512,1,2,0.75], #0.1%
                            #[512,1,2,0.69], #1%
                            [512,1,2,0.62]  #10%
                            ])
    #Epsilon is chosen by visual inspection of graphs

    results =onePercentageErrorParameters(parameters)#np.vstack([results,onePercentageErrorParameters(parameters)])

    titles = ["$n$","$\\ell$","$m$","$p$","$r$","$t$","$\\alpha$","$\\epsilon$","$c$","$a$","$\\delta$","Lattice dim. \\\\ in attack","Public key size","Blowup factor","Error prob. \\\\ upper bound"]
    #formats = ["%d","%d", "%d","%d","%d","%d","%.2f","%.2f","%.2f","%.2e","%.2f","%.2e"]
    formats = ["%8d","%8d", "%8d","%8d","%8d","%8d","%8.4f","%8.2f","%8.2f","%8.2f","%8.4f","%5d","%8.2e","%8.2f","%8.3e"]
    titles = np.array(titles)
    titlesDict = dict(zip(titles,formats))
    titlesDict.update(dict(zip(["Empirical letter \\\\ error prob.","Empirical message \\\\ error prob.","Encryption Count","Key reuse"],["%8.3e","%8.3e","%d","%d"])))
    titlesDict.update(dict(zip(["Empirical FFP-NG \\\\ error prob.","Sample Count"],["%8.3e","%d"])))

    '''
    n=256
    t=2
    r=79
    res=parameterInstation(n,t,r)
    results = np.vstack([results,res])

    n=256
    t=32
    r=1
    res=parameterInstation(n,t,r)
    results = np.vstack([results,res])
    '''
    
    elaborate=False
    generateKey=True
    
    findErrorProb = False
    findForcedErrorProb = False
    

    if(results.ndim == 1):
        results = np.reshape(results,(-1,results.size))

    #results = np.reshape(results[-1,:],(1,-1)) #Only consider the last set of parameters
    #results = np.reshape(results[0:2,:],(2,-1)) #Only consider the first two sets of parameters
    
    #results = np.reshape(results[-3:,:],(3,-1)) #Only consider the last three sets of parameters
    if findErrorProb:
        repetitions=1000 #number of keys #takes about 90 minutes with 10.000 without repeated key generation
        encryptionsPrKey = 10 #number of messages to be encrypted for each key
        expError = empiricalError(results,repetitions,encryptionsPrKey,generateKey,elaborate)
        results=np.hstack((results,np.array(expError).reshape(-1,4)))
        titles=np.append(titles, ["Empirical letter \\\\ error prob.","Empirical message \\\\ error prob.","Encryption Count","Key reuse"])
        print(titles)
        #titles.append("Empirical letter \\\\ error prob.")
        #titles.append("Empirical message \\\\ error prob.")
        #titles.append("Encryption Count")
        #titles.append("Key reuse")
        formats.append("%8.3e")
        formats.append("%8.3e")
        formats.append("%d")
        formats.append("%d")
    
    if findForcedErrorProb:
        repetitions=1000000
        err = forcedError(results,repetitions,generateKey,elaborate)
        results=np.hstack((results,np.array(err).reshape(-1,2)))
        titles = np.append(titles,["Empirical FFP-NG \\\\ error prob.","Sample Count"])
        #titles.append("Empirical FFP-NG \\\\ error prob.")
        #titles.append("Sample Count")
        formats.append("%8.3e")
        formats.append("%d")
    
    results=results.T
    #print(results)
    
    #results,titles= loadResults("results2/tmpResults220602-14,01,34.npy")
    #results,titles= loadResults("results2/combinedResults.npy")
    
    #Prepare storing results
    timestr = time.strftime("%y%m%d-%H,%M,%S")
    tmpResultsFile = "tmpResults" + timestr + ".npy" #This is a temporary file used for results from this run only
    resultFile = "combinedResults.npy" #File where results from every run is combined
    tmpLatexFile = "tmpLatex - " + timestr + ".txt" #Latex table version of result from this run only
    latexFile = "combinedLatex.txt" #Latex table version of results from every run combined
    
    folder = "results2/"
    if not os.path.exists(folder):
        os.makedirs(folder)
    tmpResultsPath = os.path.join(folder,tmpResultsFile)
    resultPath = os.path.join(folder,resultFile) #combined results
    tmpLatexPath = os.path.join(folder,tmpLatexFile)
    latexPath = os.path.join(folder,latexFile) #combiend results
    
    #Store temporary results before combinining
    saveResults(results,titles,tmpResultsPath)
    table = matrixToLatex(results,titles,titlesDict)
    print(table)
    with open(tmpLatexPath,"w") as f:
        f.write(table)
    
    #Combine results
    if(os.path.exists(resultPath)):
        results,titles = mergeResultsWithFile(results,titles,resultPath)
    #results[[11,12,13,14,15,16,17,18]]=results[[18,11,12,13,14,15,16,17]]
    #titles[[11,12,13,14,15,16,17,18]]=titles[[18,11,12,13,14,15,16,17]]
    
    
    #results[:,[4,5,6,7,8,9]]=results[:,[7,8,9,4,5,6]]
    #Store combined results
    #saveResults(results,titles,resultPath) #removed
    #saveResults(results,titles,"results2/finalTables")
    #results=results[:,4:]
    #results=results[:,[0,1,2,3,7,8,9]]
    #results=results[:,4:]
    table = matrixToLatex(results,titles,titlesDict)
    print(table)
    #with open(latexPath,"w") as f: #removed
    #    f.write(table) #removed
    
    
    #results = np.vstack([titles,results])
    
    
    #results = matrixToLatex(results,titles,formats)
    #print(results)
    
    
    #with open(filePath,"w") as f:
    #    f.write(results)
    
if __name__ == "__main__":
    main()





#np.savetxt("mydata.csv", results, delimiter=' & ', fmt='%s & %.2f & %.2f & %.2f', newline=' \\\\\n\hline\n')

"""
system=[3*(n+l)*log(p)/log(2*r+1)-m, 2*t*m*r-p]
#system=[log(p)-m, 2*t*m*r-p]
variables=[p,m]
res = nonlinsolve(system,variables)
print(res)

#system=[3*256*log(p)/log(3)-m, 2*t*m*r-p]
#res = nonlinsolve(system,variables)
#print(res)

system=[log(x,2)-p,x-p]
res = nonlinsolve(system,[x,p])
print(res)

eq1 = Eq(3*(n+l)*log(p)/log(2*r+1),m)
eq2 = Eq(2*t*m*r,p)
variables=[p,m]
res = solve([eq1,eq2],variables,dict=True)
print(res)
"""
