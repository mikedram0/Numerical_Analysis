import numpy as np
import cmath
import math

def fillv(v):
    n = len(v)
    for k in range(n//2):
        v[k] = 2.0 / (1-4.0*k*k) -1.0 / (n*n-1.0 + n%2)

    v[n//2] = (n-3.0)/( n - n%2 -1.0) \
        - 1.0 + 1.0 / (n*n-1.0 + n%2) * ((2 - n%2)*n-1.0)

    for k in range(1,(n-1)//2+1):
        v[n-k] = v[k]



def isPowerOf2(n):
    return n!=0 and not (n & (n - 1))


# with auxiliary vectors
def  fft(f, C):
    n = len(f)

    if not isPowerOf2(n):
        sys.exit(-1)
        
    if n == 1:
        C[0] = f[0]
        return
    
    fhalf = np.ndarray(n//2)
    ce = np.ndarray(n//2, dtype=complex)
    co = np.ndarray(n//2, dtype=complex)

    for i in range(n//2): fhalf[i] = f[2*i]
    fft(fhalf, ce)

    for i in range(n//2): fhalf[i] = f[2*i+1]
    fft(fhalf, co)

    pol = cmath.rect(1.0, -2.0 * math.pi / n)
    coef = complex(1.0, 0.0)

    for m in range(n//2):
        co[m] *= coef	
        C[m] = 0.5 * (ce[m] + co[m])
        C[m+n//2] = C[m] - co[m]

        coef *= pol

    return


    
def coefs(w):
    n = len(w)-1

    f = np.ndarray(n)
    fillv(f)
    fft(f,w)
    w[n] = w[0]





#  \int_{-2}^{2} \frac{1}{1+x^2} \D x = \int_{-1}^{1} \frac{2}{1+4x^2} \D x

def f(x):
    return 2.0 / (1.0+4.0*x*x)



def calc(n):
    w = np.ndarray(n+1, dtype=complex)
    coefs(w)
    
    s = complex(0.0)
    for i in range(n+1):
        x = math.cos(i*math.pi/n)
        s += w[i] * f(x)
    
    return s.real




def main():
    correct = 2.0*math.atan(2.0)
     
    n = 2
    while True:
        s = calc(n)
        diff = s-correct
        print (s, correct, diff)
        if abs(diff) < 1e-12:
            break
        n = n * 2

    print (n)



main()










