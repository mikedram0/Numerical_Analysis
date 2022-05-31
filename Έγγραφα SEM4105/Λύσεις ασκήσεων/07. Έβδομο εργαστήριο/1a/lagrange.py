import numpy

def ell(i, x, xbar):
    n = len(x) - 1
    r = 1.0

    for j in range(0, i):
        r *= (xbar-x[j]) / (x[i] - x[j])

    for j in range(i+1, n+1):
        r *= (xbar-x[j]) / (x[i] - x[j])

    return r


def p(x,y, xbar):
    n = len(x) - 1
    r = 0.0
    for  i in range(0,n+1):
        r += ell(i,x,xbar) * y[i]

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
