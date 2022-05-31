import numpy
import math

def normalize(a,b,k):
    n = len(b)
    
    for j in range(k,n):
        maxv = abs(b[j])
        for p in range(k,n):
            v = abs(a[j,p])
            if v > maxv: maxv = v
            
        for p in range(k,n): a[j,p] /= maxv       
        b[j] /= maxv
    



def swap(a,b,k,p):
    "Swap"

    n = len(b)
    
    for j in range(k,n):
        a[p,j], a[k,j] = a[k,j], a[p,j]

    b[p], b[k] = b[k], b[p]




def pivot(a,b,k):
    "Pivot"
    
    n = len(b)
    normalize(a,b,k)

    maxp = k
    for p in range(k+1,n):
        if abs(a[p,k]) > abs(a[maxp,k]):  maxp = p

    if maxp != k: swap(a,b,k,maxp)



def triang(a,b):
    "Τριγωνοποίηση"

    n = len(b)

    for k in range(n-1):
        pivot(a,b,k)

        for i in range(k+1,n):
            ell = -a[i,k]/a[k,k]

            for j in range(k,n):
                a[i,j] += ell * a[k,j]

            b[i] += ell * b[k]
        


def backsub(a,b,x):
    "Οπισθοδρόμηση"
    
    n = len(b)
    for i in range(n-1,-1,-1):
        x[i] = b[i]
        for j in range(i+1,n):
            x[i] -= a[i,j] * x[j]        
        x[i] /= a[i,i]


def gauss(a,b,x):
    triang(a,b)
    backsub(a,b,x)



def sys(x,a,b,A,B):
    n = len(x)
    for i in range(n):
        for j in range(n):
            A[i,j] = pow(x[j],i)
	
        B[i] = (pow(b,i+1)-pow(a,i+1))/(i+1)
    

        
def f(x):
    return pow(x,3) * math.sin(math.pi * x)



def main():
    x = numpy.array([-0.9,-0.7,-0.4,0.1,0.4,0.8,0.9])
    n = len(x)
    
    a = numpy.ndarray(shape=(n,n))
    b = numpy.ndarray(shape=n)
    w = numpy.ndarray(shape=n)
    
    alimit = -1.0;
    blimit = 1.0;
    
    sys(x, alimit, blimit, a, b)

    gauss(a,b,w)

    s = 0.0
    for i in range(n):
        s = s + w[i] * f(x[i])

    correct = (2.0 -12.0 /(math.pi*math.pi)) / math.pi

    print(s, correct)


main()
