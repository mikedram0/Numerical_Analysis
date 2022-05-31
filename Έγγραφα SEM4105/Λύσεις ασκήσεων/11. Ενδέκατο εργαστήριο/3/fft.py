import numpy as np
import cmath
import math
import sys

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



    

def main():
    n = 1024

    f = np.ndarray(n)

    for i in range(n):	f[i] = -0.5 + i / n

    C = np.ndarray(n, dtype=complex)

    fft(f,C)

    out = open("fourier.txt", "w")
    
    for i in range(n):
        if i==0:
            correct = 0.0
        else:
            correct = complex(0.0, 1.0/(2.0*i*math.pi))

        print(C[i], correct, file=out)


    out.close()
    
    return


main()
