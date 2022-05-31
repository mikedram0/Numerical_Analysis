import numpy as np
import math

def c(i,n):
   return 1.0 if i%n==0 else 2.0

def b(j,n):
   return 1.0 if j%(n//2)==0 else 2.0 




def coefs(w):
    n = len(w)-1

    for i in range(n+1):
        w[i] = 0.0
        for j in range(n//2+1):
            w[i] += b(j,n)  / (1-4*j*j) * math.cos(2.0*i*j*np.pi/n)
        w[i] *= c(i,n)/n



# \int_{-2}^{2} \frac{1}{1+x^2} \D x = \int_{-1}^{1} \frac{2}{1+4x^2} \D x
def f(x):
    return 2.0 / (1.0+4.0*x*x)



def calc(n,f):
    w=np.ndarray(n+1)
    coefs(w)
    s = 0.0
    for i in range(n+1):
        x = math.cos(i*np.pi/n)
        s += w[i] * f(x)

    return s


def main():

    correct = 2.0*math.atan(2.0)
    
    n = 3
    while True:
        s = calc(n,f)
        diff = s-correct
        print (s, correct, diff)
        if abs(diff) < 1e-12: break
        n=n+1

    print (n)



main()
