'''
By Andreas H. Jørgensen (s202743)

THE LWE Encryption system
This version uses 64 bit integers to avoid overflow problems
'''

import numpy as np
import time
np.set_printoptions(suppress=True)

class LWE:
    def __init__(self,n,m,l,t,r,p,alpha,elaborate=False):
        '''
        Intializing object for encryption
        '''
        self.n=n
        self.m=m
        self.l=l
        self.t=t
        self.r=r
        self.p=p
        self.alpha=alpha
        self.elaborate=elaborate
    
    def gen(self,keyname=None,saveKey=True):
        '''
        This corresponds to the GEN-function of the system
        '''
        if(self.elaborate):
            print("\n##########\nGenerating\n########")
        #Private key
        S = np.random.randint(0,self.p,(self.n,self.l),dtype=np.int64)
        Kpriv = S
        if(self.elaborate):
            print("S", S.shape)
            print(S)
            print(S.dtype)
        
        #Public key
        A = np.random.randint(0,self.p,(self.n,self.m),dtype=np.int64)
        std = self.alpha*self.p/np.sqrt(2*np.pi)
        E = np.mod(np.random.normal(0,std,(self.l,self.m)).round(),self.p)
        B=np.mod(np.matmul(S.T,A)+E,self.p)

        Kpub = (A,B)
        if self.elaborate:
            print("A", A.shape)
            print(A)
            print("E", E.shape)
            print(E)
            print("B", B.shape)
            print(B)
        if keyname==None:
            keyname ="key"+str(self.m)+".npy"
        if saveKey:
            with open(keyname,"wb") as f:
                np.save(f,Kpriv)
                np.save(f,Kpub[0])
                np.save(f,Kpub[1])
        return (Kpriv,Kpub)
    
    def enc(self,v,Kpub):
        '''
        Function for doing encryption
        '''
        if self.elaborate:
            print("\n##########\nEncrypting\n########")
        A,B = Kpub
        n,m = A.shape
        l,m2 = B.shape
        l2, = v.shape
        if m!=self.m or m2!=self.m or n!=self.n or l!=self.l or l2!=self.l:
            print("Some error in shape of input during encryption")
            return None
        a = np.random.randint(-self.r,self.r+1,m)
        if self.elaborate:
            print("a", a.shape)
            print(a)
        u = np.mod(np.matmul(A,a),self.p)
        c = np.mod(np.matmul(B,a)+self.f(v),self.p)
        if(self.elaborate):
            print(c.shape)
        return (u,c)
    
    def dec(self,cipher,Kpriv):
        '''
        Function for doing decryption
        '''
        if self.elaborate:
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
        return self.finv(np.mod(c-np.matmul(S.T,u),self.p))
        
    def f(self,v):
        return (self.p/self.t*v).round()

    def finv(self,c):
        return (self.t/self.p*c).round()
    

def main():
    elaborate = False
    generateKey = True
    ##Lattice based cryptography article
    #n=233; l=233; m=4536; p=32749; r=1; t=40; alpha = 0.000217
    #n=136; l=136; m=2008; p=2003;  r=1; t=2; alpha = 0.0065
    #n=166; l=166; m=1319; p=4093; r=4; t=2; alpha = 0.0024
    #n=192; l=192; m=1500; p=8191; r=5; t=4; alpha = 0.0009959
    #My own parameters (See word)
    #n=136; l=136; m=8000;  p=34000; r=1; t=2; alpha=0.006321
    #n=136; l=136; m=8001;  p=34000; r=1; t=2; alpha=0.0053
    #n=136; l=136; m=8002;  p=34000; r=1; t=2; alpha=0.001528265
    n=128; l=128; m=7176; p=28703; r=1; t=2; alpha=0.0008
    
    np.random.seed(seed=42) #control randomness
    
    ##Preparing/setup
    lwe = LWE(n,m,l,t,r,p,alpha,elaborate)
    if generateKey:
        keys = lwe.gen()
        priv,pub = keys
    else:
        with open("key"+str(m)+".npy","rb") as f:
            priv = np.load(f)
            pub = np.load(f),np.load(f)
    
    #Encryption
    v = np.random.randint(0,t,l)
    print("Message:", v)
    cipher = lwe.enc(v,pub)
    if cipher == None:
        return
    print("Ciphertext:",cipher)
    
    #Decryption
    m = lwe.dec(cipher,priv)
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