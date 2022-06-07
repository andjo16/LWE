'''
THE LWE Encryption system
This is the version that is working
'''

import numpy as np
import time
np.set_printoptions(suppress=True)

class LWE:
    def __init__(self,n,m,l,t,r,q,alpha):
        '''
        Intializing object for encryption
        '''
        self.n=n
        self.m=m
        self.l=l
        self.t=t
        self.r=r
        self.q=q
        self.alpha=alpha
    
    def gen(self,keyname=None,elaborate=False):
        '''
        This corresponds to the GEN-function of the system
        '''
        print("\n##########\nGenerating\n########")
        #Private key
        S = np.random.randint(0,self.q,(self.n,self.l))
        #S=np.array([[5]*self.l]*self.n).T
        Kpriv = S
        print("S", S.shape)
        print(S)
        print(S.dtype)
        
        #Public key
        A = np.random.randint(0,self.q,(self.m,self.n))
        std = self.alpha*self.q/np.sqrt(2*np.pi)
        E = np.mod(np.random.normal(0,std,(self.m,self.l)).round(),self.q)
        #print("E", E.shape)
        #print(E)
        P=np.mod(matmulmod(A,S,self.q)+E,self.q)
        print(np.mod(dotmod(A[0,:],S[:,0],self.q)+E[0,0],self.q))

        Kpub = (A,P)
        if elaborate:
            print("A", A.shape)
            print(A)
            print("E", E.shape)
            print(E)
            print("P", P.shape)
            print(P)
        if keyname==None:
            keyname ="key"+str(self.m)+".npy"
        with open(keyname,"wb") as f:
            np.save(f,Kpriv)
            np.save(f,Kpub[0])
            np.save(f,Kpub[1])
        return (Kpriv,Kpub)
    
    def enc(self,v,Kpub, elaborate=False):
        '''
        Function for doing encryption
        '''
        print("\n##########\nEncrypting\n########")
        A,P = Kpub
        #print(A.shape)
        m,n = A.shape
        m2,l = P.shape
        l2, = v.shape
        if m!=self.m or m2!=self.m or n!=self.n or l!=self.l or l2!=self.l:
            print("Some error in shape of input during encryption")
            return None
        a = np.random.randint(-self.r,self.r+1,m)
        if elaborate:
            print("a", a.shape)
            print(a)
        u = np.mod(matmulmod(A.T,a,self.q),self.q)
        c = np.mod(matmulmod(P.T,a,self.q)+self.f(v),self.q)
        print(c.shape)
        return (u,c)
    
    def dec(self,cipher,Kpriv,elaborate=False):
        '''
        Function for doing decryption
        '''
        print("\n##########\nDecrypting\n########")
        u,c = cipher
        S = Kpriv
        n,l = S.shape
        n2, = u.shape
        l2, = c.shape
        if n!=self.n or n2!=self.n or l!=self.l or l2 !=self.l:
            print("Some error in shape of input during Decryption")
            print("%d: %d, %d || %d: %d, %d" %(self.n,n,n2,self.l,l,l2))
            return None
        return self.finv(np.mod(c-matmulmod(S.T,u,self.q),self.q))
        
    def f(self,v):
        print(v,self.q,self.t)
        print(type(v),v.dtype,type(self.q),type(self.t))
        return (int(self.q)/int(self.t)*v).round()

    def finv(self,c):
        print(c,self.q,self.t)
        print(type(c),c.dtype,type(self.q),type(self.t))
        return (self.t/self.q*c).round()


matmulmodTime = 0
def dotmod(u,v,q,x=-1,y=-1):
    #print(u)
    #print(v)
    if y==0 and x%50==0:
        tmp = (x,y)
        print(tmp, end="; ")
        print("time Elapsed: %.2f" %(time.time()-matmulmodTime))
    res = 0
    for i in range(len(u)):
        res += np.mod(u[i]*v[i],q)
    return res

def matmulmod(A,B,q):
    global matmulmodTime
    matmulmodTime = time.time()
    if(A.ndim>2):
        print("First matrix has %d>2 dimensions" % A.ndim)
        raise ValueError()
    elif(A.ndim==1):
        k1,=A.shape
        A = A.reshape((1,k1))
    if B.ndim>2:
        print("First matrix has %d>2 dimensions" % B.ndim)
        raise ValueError()
    if(B.ndim==1):
        print("HERE: Second array has one dimension")
        k2,=B.shape
        B = B.reshape((k2,1))
    m,k1 = A.shape
    k2,n = B.shape
    if(k1!=k2):
        print("%d columns of first matrix does not match %d rows of second matrix" %(k1,k2))
        raise ValueError()
    return np.array([[dotmod(A[i,:],B[:,j],q,i,j) for j in range(n)] for i in range(m)]).squeeze()
        

def main():
    elaborate = True
    generateKey = True
    ##Lattice based cryptography article
    #n=233; l=233; m=4536; q=32749; r=1; t=40; alpha = 0.000217
    n=136; l=136; m=2008; q=2003;  r=1; t=2; alpha = 0.0065
    #n=166; l=166; m=1319; q=4093; r=4; t=2; alpha = 0.0024
    #n=192; l=192; m=1500; q=8191; r=5; t=4; alpha = 0.0009959
    
    #My own parameters (See word)
    #n=136; l=136; m=8000;  q=34000; r=1; t=2; alpha=0.006321
    #n=136; l=136; m=8001;  q=34000; r=1; t=2; alpha=0.0053
    #n=136; l=136; m=8002;  q=34000; r=1; t=2; alpha=0.001528265
    
    np.random.seed(seed=42) #control randomness
    
    ##Preparing/setup
    lwe = LWE(n,m,l,t,r,q,alpha)
    if generateKey:
        keys = lwe.gen(elaborate)
        priv,pub = keys
    else:
        with open("key"+str(m)+".npy","rb") as f:
            priv = np.load(f)
            pub = np.load(f),np.load(f)
    
    #Encryption
    v = np.random.randint(0,t,l)
    print("Message:", v)
    cipher = lwe.enc(v,pub,elaborate)
    if cipher == None:
        return
    print("Ciphertext:",cipher)
    
    #Decryption
    m = lwe.dec(cipher,priv,elaborate)
    if m is None:
        return
    
    #Verification of results
    print("Original message:")
    print(v)
    print("Decrypted message:")
    print(m)
    print("Decrypted message modulo t:")
    mred = np.mod(m,lwe.t)
    print(mred)
    
    equal = v==mred
    errors = np.where(equal==False)
    print("Found %d errors" %len(errors[0]))
    if(len(errors[0])<20):
        for i in errors[0]:
            print("Error at index %d: %d != %d" %(i,v[i],mred[i]))
    
    #print(v==np.mod(m,t))


if __name__ == "__main__":
    main()