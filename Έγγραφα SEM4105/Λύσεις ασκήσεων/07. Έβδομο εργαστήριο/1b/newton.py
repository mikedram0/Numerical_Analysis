import numpy



def q(i, x, xbar):
    if (i == 0):
        return 1.0
    else:
        return q(i-1,x,xbar) * (xbar - x[i-1])
    


def coef(x, y):
    n = len(x) - 1

    a = numpy.ndarray(shape=(n+1), dtype=numpy.float64)
    a[0] = y[0]

    for j in range(1,n+1):
        a[j] = y[j]
        for i in range(0,j):
            a[j] -= a[i] * q(i,x,x[j])
        a[j] /= q(j,x, x[j])

    return a


def p(x, y, xbar):
    a = coef(x,y)
    n = len(x) - 1
    r = 0.0
    for  i in range(0,n+1):
        r += q(i,x,xbar) * a[i]

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
    
    xmin = min(x)
    xmax = max(x)
    m = 100 - 1 
    h = (xmax-xmin) / m
    for i in range (0,m+1):
        xbar = xmin + i * h
        print (xbar, p(x,y,xbar), numpy.sin(xbar))

    return


        
main()
