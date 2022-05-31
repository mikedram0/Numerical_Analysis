import numpy
import math


"""
  y' = z
  z' = (x*x-5) * y   
"""
def f(x,y):
    return numpy.array([y[1], (x*x-5.0)*y[0]])


def ralston(t0, y0, t1, f):
    h = t1-t0

    k1 = h*f(t0, y0)
    
    p = y0 + 2.0/3.0 * k1
    
    k2 = h*f(t0+2.0/3.0*h, p)

    y1 = y0 + (k1 + 3.0 * k2)/4.0

    return y1



def main():
    nsteps = 101
    a = -2.0
    b = 2.0
    h = (b-a)/(nsteps-1)    
    offset = nsteps//2


    y = numpy.zeros(nsteps)
    v = numpy.zeros(nsteps)
    # nsteps/2 , 0.0, nsteps/2

    y0 = numpy.array([-1.0/math.sqrt(2.0*math.sqrt(math.pi)), 0.0])

    y[offset], v[offset] = y0[0], y0[1]

    x = numpy.zeros(nsteps)
    for j in range(nsteps):
        x[j] = a+h*j

    for j in range(nsteps//2):
        x0 = x[offset+j]
        y1 = ralston(x0, y0, x0+h, f)

        y0 = y1
        y[offset+(j+1)], v[offset+(j+1)] = y1[0], y1[1]
  

    y0 = numpy.array([-1.0/math.sqrt(2.0*math.sqrt(math.pi)), 0.0])

    for j in range(nsteps//2):
        x0 = x[offset-j]
        y1 = ralston(x0, y0, x0-h, f)
        
        y0 = y1
        
        y[offset-(j+1)], v[offset-(j+1)] = y1[0], y1[1]
  


    out = open("psi.txt", "w")
    for j in range(nsteps):
        s = str(x[j])+ ' ' + str(y[j]) + ' ' + str(v[j]) + '\n'
        out.write(s)
    
    out.close()



main()
