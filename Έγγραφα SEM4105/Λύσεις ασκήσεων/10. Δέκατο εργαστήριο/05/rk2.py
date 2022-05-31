import math


def correct(x):
    return 1.0 - math.exp(-x) + x*x -x


def f(x,y):
    return x*x + x-y

def heun(x0,y0,x1, f):
    h = x1-x0
    
    k1 = h * f(x0,y0)
    k2 = h * f(x1, y0 + k1)

    return y0 + (k1 + k2)/2.0


def ralston(x0, y0, x1, f):
    h = x1-x0
    
    k1 = h * f(x0,y0)
    k2 = h * f(x0 + 2.0/3.0 * h, y0 + 2.0 / 3.0 * k1)

    return y0 + (k1 + 3.0 * k2) / 4.0


def main():
    a = 0.0
    b = 0.6
    h = 0.01
    
    yinit = 0.0

    x = a
    yh = yinit
    yr = yinit

    
    while x < b:
        if x+h > b:
	    h = b-x
	

        yh = heun(x, yh, x+h, f)
	yr = ralston(x, yr, x+h, f)
	
	x += h

        print (x, yh, yr, correct(x));


main()
