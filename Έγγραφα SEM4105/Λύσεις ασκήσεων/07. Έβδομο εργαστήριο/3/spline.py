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
    triang(a,b)
    backsub(a,b,x)




    
"""

Για i=0,...,n-1 :   
    dx_i = x_{i+1}-x_i

Για i=0,...,n-1 :

    a_i dx_i^2 + b_i dx_i + c_i = (f_{i+1} - f_i) / dx_i


Για i=0,...,n-2 :

    3 a_i dx_i^2 + 2b_i dx_i + c_i - c_{i+1} = 0


Για i=0,...,n-2 :

    3a_i dx_i + b_i - b_{i+1} = 0



b_0 = 0

3 a_{n-1} dx_{n-1} + b_{n-1} = 0

"""

def spline_calculate(x, y, a, b, c, d):

    n = len(x)-1
    m = 3*n
    
    A = numpy.zeros(shape=(m,m), dtype=numpy.float64)
    B = numpy.zeros(shape=m, dtype=numpy.float64)

    # X = a_0, a_1, ..., a_{n-1}, b_0, b_1, ..., b_{n-1}, c_0, c_1, ..., c_{n-1}


    # Πρώτο σετ από n εξισώσεις: I από 0 έως και n-1
    I = 0
    for i in range(n):
        dx = x[i+1]-x[i]
        A[I,i] = dx*dx        # coefficient of a_i 
        A[I,i+n] = dx       # coefficient of b_i 
        A[I,i+2*n] = 1.0    # coefficient of c_i 
        B[I] = (y[i+1]-y[i]) / dx
        I=I+1

        

    # Δεύτερο σετ από  n-1 εξισώσεις: I από n έως και 2*n-2
    for i in range (n-1):
        dx = x[i+1]-x[i]

        A[I,i] = 3*dx*dx       # coefficient of a_i 
        A[I,I+n] = 2.0*dx      # coefficient of b_i 
        A[I,i+2*n] = 1.0       # coefficient of c_i 
        A[I,i+2*n+1] = -1.0   # coefficient of c_{i+1} 

        I=I+1
        

    # Τρίτο σετ από  n-1 εξισώσεις: I από 2*n-1 έως και 3*n-3
    for i in range (n-1):
        dx = x[i+1]-x[i]
        
        A[I,i] = 3.0*dx    # coefficient of a_i 
        A[I,i+n] = 1.0     # coefficient of b_i 
        A[I,i+n+1] = -1.0  # coefficient of b_{i+1}

        I=I+1
        

    # Οριακές συνθήκες:
        
    # Ι = m-2
    A[I,n] = 1.0   # coefficient of b_0
    I=I+1

    # I = m-1
    A[I,n-1] = 3.0 * (x[n]-x[n-1])   # coefficient of a_{n-1}
    A[I,2*n-1] = 1.0                 # coefficient of b_{n-1}


    X= numpy.ndarray(shape=m, dtype=numpy.float64)

    gauss(A,B,X)
    
    for i in range(n):
        a[i] = X[i]
        b[i] = X[i+n]
        c[i] = X[i+2*n]
        d[i] = y[i]
        

    return



def spline_evaluate(a,b,c,d,x,z):
    n = len(x)-1
    
    i  = 0    
    while i < n and x[i+1] < z :  i=i+1
    

    if z > x[n] : i = n-1


    dx = z-x[i]
    
    return ((a[i]*dx+b[i])*dx+c[i])*dx+d[i]



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

    n = len(x)-1

    a = numpy.ndarray(shape=n, dtype=numpy.float64)
    b = numpy.ndarray(shape=n, dtype=numpy.float64)
    c = numpy.ndarray(shape=n, dtype=numpy.float64)
    d = numpy.ndarray(shape=n, dtype=numpy.float64)

    spline_calculate(x,y,a,b,c,d)

    xmin = min(x)
    xmax = max(x)
    m = 100 - 1 
    h = (xmax-xmin) / m
    for i in range (0,m+1):
        xbar = xmin + i * h
        print(xbar, spline_evaluate(a,b,c,d,x,xbar), numpy.sin(xbar))

    return


main()
