import numpy
import math

def runge_kutta(A,b,c,x0,y0,x1, f):
    s = len(c)
    
    h = x1-x0

    k = numpy.zeros(s)

    p = y0
    
    for i in range(s):
        sum = y0
        for j in range(i):
            sum += A[i,j] * k[j]

        k[i] = h * f(x0+c[i] * h, sum)
        p += b[i] * k[i]

    return p
        

def correct(x):
    return x + math.sqrt(1+2.0*x*x)


def f(x,y):
    return (y+x)/(y-x)


def main():
    a = 0.0
    b = 2.0
    h = 0.1
    yinit = 1.0

    x = a
    y = yinit


    rk_A=numpy.array([[0.0, 0.0, 0.0, 0.0,],
	              [0.5, 0.0, 0.0, 0.0,],
	              [0.0, 0.5, 0.0, 0.0,],
	              [0.0, 0.0, 1.0, 0.0]])
                     
    rk_b =  numpy.array([1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0])

    rk_c =  numpy.array([0.0, 0.5, 0.5, 1.0])

    while x < b:
        if x+h > b:
            h = b-x
	
        y = runge_kutta(rk_A, rk_b, rk_c, x, y, x+h, f)
	
        x += h
        print (x,y,correct(x))


main()


