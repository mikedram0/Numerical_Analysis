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



def gauss(a,b,x):
    "Gauss"

    triang(a,b)
    backsub(a,b,x)




def nonlinear(x, f, a):
    f[0] = x[0]             + x[1]             + x[2] - 3.0
    f[1] = x[0]*x[0] * x[1] + x[1]*x[1] * x[2] + x[2]*x[2] * x[0] - 4.0
    f[2] = x[0]*x[0]        + x[1]*x[1]        + x[2]*x[2]  - 5.0


    a[0,0] = 1.0
    a[0,1] = 1.0
    a[0,2] = 1.0

    a[1,0] = 2.0 * x[0] * x[1] + x[2]*x[2]
    a[1,1] = 2.0 * x[1] * x[2] + x[0]*x[0]
    a[1,2] = 2.0 * x[2] * x[0] + x[1]*x[1]

    a[2,0] = 2.0 * x[0]
    a[2,1] = 2.0 * x[1]
    a[2,2] = 2.0 * x[2]

    

def check(x, tol):
    return numpy.all(abs(x) < numpy.ones(len(x))*tol)




def main():
    n = 3
    tol = 1e-7

    x = numpy.ndarray(shape=(n), dtype=numpy.float64)
    x[0] = 1.0
    x[1] = 1.5
    x[2] = 2.1
    
    y = numpy.ndarray(shape=(n), dtype=numpy.float64)
    b = numpy.ndarray(shape=(n), dtype=numpy.float64)
    a = numpy.ndarray(shape=(n,n), dtype=numpy.float64)

    while True:
        nonlinear(x,b,a)
        
        if check(b,tol): break
            
        gauss(a,b,y)
        
        x -= y

    print (x)




    
main()
