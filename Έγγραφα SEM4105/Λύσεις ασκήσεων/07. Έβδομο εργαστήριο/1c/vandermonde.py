import numpy


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


        

def pol(a, xbar):
    r = 0.0

    n = len(a)
    for i in range(0,n):
        r += a[i] * xbar**i
    
    return r



        

def readxy(filename):
    f = open(filename, "r")
    n = int(f.readline())
    
    x = numpy.ndarray(shape=n, dtype=numpy.float64)
    y = numpy.ndarray(shape=n, dtype=numpy.float64)

    for i in range(n):
        x[i], y[i] = numpy.fromfile(f, dtype=numpy.float64, count=2, sep=" ")

    return x,y



def main():

    x,y = readxy("points.dat")
    n = len(x)

    a = numpy.ndarray(shape=(n,n), dtype=numpy.float64)

    for i in range (0,n):
        for j in range (0,n):
            a[i,j] = x[i]**j
    
    coef = numpy.ndarray(shape=(n), dtype=numpy.float64)

    triang(a,y)
    backsub(a,y,coef)


    xmin = min(x)
    xmax = max(x)
    m = 100 - 1 
    h = (xmax-xmin) / m
    for i in range (0,m+1):
        xbar = xmin + i * h
        print (xbar, pol(coef,xbar), numpy.sin(xbar))

    return


        
main()
