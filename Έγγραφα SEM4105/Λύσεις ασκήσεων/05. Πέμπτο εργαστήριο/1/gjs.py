import numpy


def check(x, tol):
    return numpy.all(abs(x) < numpy.ones(len(x))*tol)

def gauss_jacobi(a,b,x):
    n = len(b)

    x.fill(0.0)
    
    while True:
        y = b - numpy.matmul(a, x)

        for i in range(0,n):
            y[i] /= a[i,i]

        x += y

        if check(y,1e-7):  break
    


def gauss_seidel(a,b,x):
    n = len(b)

    x.fill(0.0)
    
    while True:
        y = b - numpy.matmul(a, x)

        for i in range(0,n):
            y[i] /= a[i,i]
            x[i] += y[i]

        if check(y,1e-7):  break
    


def main():
    n = 4
    tol = 1e-7

    a = numpy.ndarray(shape=(n,n), dtype=numpy.float64)
    b = numpy.ndarray(shape=(n), dtype=numpy.float64)
    x = numpy.ndarray(shape=(n), dtype=numpy.float64)

    a[0,0] = 12.1
    a[0,1] = 3.9
    a[0,2] = 0.3
    a[0,3] = -4.1

    a[1,0] = 4.3
    a[1,1] = -11.3
    a[1,2] = 0.8
    a[1,3] = 1.5

    a[2,0] = 1.0
    a[2,1] = -2.8
    a[2,2] = 14.3
    a[2,3] = -8.1

    a[3,0] = 2.4
    a[3,1] = 6.1
    a[3,2] = -1.1
    a[3,3] = 12.5

    b[0] = 1.2
    b[1] = 2.3
    b[2] = 3.4
    b[3] = 4.5


    gauss_jacobi(a,b,x)
    
    print ("solution with jacobi:\n")
    for i in range(0,n):  print (x[i])
    print ("\n")

    gauss_seidel(a,b,x)
    
    print ("solution with seidel:\n")
    for i in range(0,n):  print (x[i])
    print ("\n")


    
    
main()
