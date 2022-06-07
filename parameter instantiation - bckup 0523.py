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
    Slightly alternative way to find p and m
    '''
    #First we find p and m based on n,l, t and the lower bounds on p and m
    p,m=symbols('p, m', real=True)
    #eq1 = Eq(3*(n+l)*log(p)/log(2*r+1),m) #remember to update m
    #eq1 = Eq((2*(n+l)*log(p)+200)/log(2*r+1)-m)
    #eq1 = Eq((2*(n+l)+200)*log(p)/log(2*r+1),m)
    #eq2 = Eq(2*t*m*r-p)
    #eqs = [(2*(n+l)*log(p)+200)/log(2*r+1)-m,2*t*m*r-p]
    
    def f(p):
        return (2*(n+l)*math.log(p[0],2)+200)/math.log(2*r+1,2)-p[0]/(2*t*r)
    
    p=fsolve(f,[10000])
    print(p)
    if(p>10):
            #p needs to be a prime
            p=math.ceil(p)
            #m=math.ceil(sol[m])
            if not sympy.isprime(p):
                #print(type(p))
                p=sympy.nextprime(p)
                
                #p=2*t*m*r
            p=int(p)
            #Update m according to p
            #m=math.ceil(3*(n+l)*math.log(p)/math.log(2*r+1)) #Remember to update m
            m=math.ceil((2*(n+l)*math.log(p,2)+200)/math.log(2*r+1,2))
    else:
        print("WARNING: p is very small")
    '''
    #eqs = (2*(n+l)*log(p,2))/log(2*r+1,2)-p/(2*t*r)
    #variables=[p]
    #res = solve(eqs,variables,dict=True)
    #print(res)
    #res = [{x: si[x].evalf() for x in si } for si in res]
    #print(res)
    for sol in res:
        #print(type(sol))
        #print(sol.keys())

        if(sol[p]>10):
            #p needs to be a prime
            p=math.ceil(sol[p])
            #m=math.ceil(sol[m])
            if not sympy.isprime(p):
                print(type(p))
                p=sympy.nextprime(p)
                
                #p=2*t*m*r
            p=int(p)
            #Update m according to p
            #m=math.ceil(3*(n+l)*math.log(p)/math.log(2*r+1)) #Remember to update m
            m=math.ceil((2*(n+l)*math.log(p)+200)/math.log(2*r+1))
            break
    '''
    return p,m


def findPandMold(n,l,t,r):
    #First we find p and m based on n,l, t and the lower bounds on p and m
    p,m=symbols('p, m', real=True)
    #eq1 = Eq(3*(n+l)*log(p)/log(2*r+1),m) #remember to update m
    eq1 = Eq((2*(n+l)*log(p)+200)/log(2*r+1),m)
    #eq1 = Eq((2*(n+l)+200)*log(p)/log(2*r+1),m)
    eq2 = Eq(2*t*m*r,p)
    variables=[p,m]
    res = solve([eq1,eq2],variables, dict=True, manual=True)
    res = [{x: si[x].evalf() for x in si } for si in res]
    #print(res)
    for sol in res:
        #print(type(sol))
        #print(sol.keys())

        if(sol[p]>10):
            #p needs to be a prime
            p=math.ceil(sol[p])
            #m=math.ceil(sol[m])
            if not sympy.isprime(p):
                print(type(p))
                p=sympy.nextprime(p)
                
                #p=2*t*m*r
            p=int(p)
            #Update m according to p
            #m=math.ceil(3*(n+l)*math.log(p)/math.log(2*r+1)) #Remember to update m
            m=math.ceil((2*(n+l)*math.log(p)+200)/math.log(2*r+1))
            break
    return p,m


def parameterInstation(n=128,t=2,r=1):
    l=n
    p,m=findPandM(n,l,t,r)

    #Now we find c (the exponent of gamma, as p>zeta*sqrt(ln(n)/n)
    c = math.log(p/(2*math.sqrt(math.log(n)/n)))/math.log(n)
    c=math.floor(c * 100)/100.0


    #epsilon = Symbol("epsilon", real=True)
    #system = [n/((n**c)*math.sqrt(math.log(n)))-1/(t*math.sqrt((r*(r+1)*m)/(6*math.pi))*(math.log(n))**(1/2+epsilon))]
    #variables = [epsilon]
    #res = nonlinsolve(system,variables)
    #print(res)
    
    #We can now find epsilon as the largest value that allows the upper and lower bounds on alpha to meet
    epsilon = math.log(n**c*math.sqrt(math.log(n))/(n*t*math.sqrt((r*(r+1)*m)/(6*math.pi))))/(math.log(math.log(n)))-1/2
    epsilon=math.floor(epsilon * 100)/100.0
    #print(epsilon)
    
    alpha = 1/(t*math.sqrt((r*(r+1)*m/(6*math.pi)))*(math.log(n))**(1/2+epsilon))
    
    #Find the value of a that minimizes error
    def f(a):
        return l*2*math.exp(-2*m*r**2/a[0]**2)+l/((a[0]-1)/a[0]*(math.log(n)**(1/2+epsilon))/2)*math.exp(-(((a[0]-1)/a[0]*(math.log(n))**(1/2+epsilon))**2)/8)
        #return l*2*math.exp(-2*m*r^2/a[0]**2)#+l/((a-1)/a[0]*(math.log(n)**(1/2+epsilon)))*math.exp(-(((a[0]-1)/a[0])**2*(math.log(n))**(1/2+epsilon))/8)
    '''
    a = Symbol('a')
    function = l*2*exp(-2*m*r/a**2)+l/((a-1)/a*(sympy.log(n)**(1/2+epsilon)))*sympy.exp(-(((a-1)/a)**2*(sympy.log(n))**(1/2+epsilon))/8)
    #function = l/((a-1)/a*(sym.log(n)**(1/2+epsilon)))*exp((((a-1)/a)**2*(sym.log(n))**(1/2+epsilon))/8)
    #function = l/((a-1)/a*(sympy.log(n)**(1/2+epsilon)))*sympy.exp((((a-1)/a)**2*(sympy.log(n))**(1/2+epsilon))/8)
    #mini = minimum(function,a)
    #print(mini)
    
    
    lower_bound = 2
    upper_bound = 20
    zeros = solveset(function, a, domain=Interval(lower_bound, upper_bound))
    print(zeros)
    assert zeros.is_FiniteSet # If there are infinite solutions the next line will hang.
    ans = Min(function.subs(a, lower_bound), function.subs(a, upper_bound), *[function.subs(a, i) for i in zeros])
    
    print(ans)
    '''
    
    a = fmin(f,np.array([10]),disp=False)
    
    
    keySize = math.ceil((n+l)*m*math.log(p))
    ciphertextBlowup = (1+n/l)*math.log(p)/math.log(t)
    #print("n,l,m,t,r,p,epsilon,c,a,key,blowup,error")
    #print(n,l,m,t,r,p,epsilon,c,a[0],keySize,ciphertextBlowup, f(a))
    
    #Return all results 
    #print(f(a),f(a)*100)
    ress = np.array([n,l,m,p,r,t,alpha,epsilon,c,a[0],keySize,ciphertextBlowup,f(a)*100])
    print("types",type(t),type(p))
    return ress

def onePercentageErrorParameters():
                             #n,r,t,epsilon
    parameters = np.array([[256,1,2,0.82],
                            [256,1,2,0.75],
                            [256,1,2,0.67]])
    #Epsilon is chosen by visual inspection of graphs
    ress = np.empty((0,13))
    for params in parameters:
        n,r,t,epsilon = params
        l=n
    
        p,m = findPandM(n,l,t,r)
        
        #Now we find c (the exponent of gamma, as p>zeta*sqrt(ln(n)/n)
        c = math.log(p/(2*math.sqrt(math.log(n)/n)))/math.log(n)
        c=math.floor(c * 100)/100.0
        
        alpha = 1/(t*math.sqrt((r*(r+1)*m/(6*math.pi)))*(math.log(n))**(1/2+epsilon))
        keySize = math.ceil((n+l)*m*math.log(p))
        ciphertextBlowup = (1+n/l)*math.log(p)/math.log(t)
        
        def f(a):
            return l*2*math.exp(-2*m*r**2/a[0]**2)+l/((a[0]-1)/a[0]*(math.log(n)**(1/2+epsilon))/2)*math.exp(-(((a[0]-1)/a[0]*(math.log(n))**(1/2+epsilon))**2)/8)
        a = fmin(f,np.array([10]),disp=False)
        print(f(a))
        
        ress = np.vstack([ress,[n,l,m,p,r,t,alpha,epsilon,c,a[0],keySize,ciphertextBlowup,f(a)*100]])
    return ress


def testDecryptErrorRateProcess(n,l,m,p,r,t,alpha,encCount=1,elaborate=False):
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

def testDecryptionErrorRateMultiprocessing(n,l,m,p,r,t,alpha,repetitions=100,encCount=1,generateKey=False,elaborate=False):
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
            if (i+1)%50==0:
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

def testDecryptionErrorRate(n,l,m,p,r,t,alpha,repetitions=100,generateKey=False,elaborate=False):
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
        experimentalError.append([letterFraction,msgFraction])
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
        experimentalError.append(FFPNG.errorProbParallel(m,p,t,r,alpha,repetitions,elaborate))
        #print("time: %0.2f" %(time.time()-t))
    return experimentalError


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
        print(i)
        if i<=9:
            tmp = "& %15s " %titles[i]#%row[0]
        else:
            tmp = "\\multicolumn{2}{|c|}{\\makecell{%15s}}" %titles[i]
        for d in row:
            tmp += " & "+formats[i] %d
            if i in [12,13,14,15]:
                tmp += "\\%"
        if i ==6 or i>=9:
            tmp += "\\\\\n\hline\n"
            if i==6:
                tmp +="\multirow{3}{*}{\\shortstack{Meta\\\\ parameters}}\n"
            #if i==9:
            #    tmp +="\multirow{3}{*}{}\n"
        else:
            tmp+= "\\\\\n\cline{2-%d}\n" %(columns+2)
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
    #results = np.array(["$n$","$\ell$","$m$","$p$","$r$","$t$","$\\alpha$","$\epsilon$","$c$","$a$","Public key size","Blowup factor","Error prob."])
    n=128
    t=2
    r=1
    res = parameterInstation(n,t,r)
    results =res #np.vstack([results,res])
    n=256
    res=parameterInstation(n,t,r)
    results = np.vstack([results,res])
    n=512
    res=parameterInstation(n,t,r)
    results = np.vstack([results,res])
    
    titles = ["$n$","$\ell$","$m$","$p$","$r$","$t$","$\\alpha$","$\epsilon$","$c$","$a$","Public key size","Blowup factor","Error prob. \\\\ upper bound"]
    #formats = ["%d","%d", "%d","%d","%d","%d","%.2f","%.2f","%.2f","%.2e","%.2f","%.2e"]
    formats = ["%8d","%8d", "%8d","%8d","%8d","%8d","%8.4f","%8.2f","%8.2f","%8.2f","%8.2e","%8.2f","%8.3e"]
    
    results = np.vstack([results,onePercentageErrorParameters()])
    
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
    
   
    #print(results.shape)
    #print(results[-1,:])
    #results = np.reshape(results[-1,:],(1,-1)) #Only consider the last set of parameters
    #results = np.reshape(results[0:2,:],(2,-1)) #Only consider the first two sets of parameters
    
    #results = np.reshape(results[-3:,:],(3,-1)) #Only consider the last three sets of parameters
    if findErrorProb:
        repetitions=1000 #number of keys #takes about 90 minutes with 10.000 without repeated key generation
        encryptionsPrKey = 1 #number of messages to be encrypted for each key
        expError = empiricalError(results,repetitions,encryptionsPrKey,generateKey,elaborate)
        results=np.hstack((results,np.array(expError).reshape(-1,2)))
        titles.append("Empirical letter \\\\ error prob.")
        titles.append("Empirical message \\\\ error prob.")
        formats.append("%8.3e")
        formats.append("%8.3e")
    
    if findForcedErrorProb:
        repetitions=10000
        err = forcedError(results,repetitions,generateKey,elaborate)
        results=np.hstack((results,np.array(err).reshape(-1,1)))
        titles.append("Empiical FFP-NG \\\\ error prob.")
        formats.append("%8.3e")
    
    #results = np.vstack([titles,results])
    results=results.T
    print(results)
    results = matrixToLatex(results,titles,formats)
    print(results)
    
    timestr = time.strftime("%y%m%d-%H,%M,%S")
    fileName = "LatexTable - " + timestr + ".txt"
    filePath = os.path.join('results', fileName)
    with open(filePath,"w") as f:
        f.write(results)
    
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
