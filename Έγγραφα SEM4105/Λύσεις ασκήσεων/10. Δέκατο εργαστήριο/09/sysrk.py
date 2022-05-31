import numpy


def force(x):
    return x * (1.0 - 0.01*x*x)


"""
  m x'' = F(x) =>    x' = v     =>    y0 = x, y1 = v
                     v' = F/m
"""
def f(t,y):
    mass = 2.0

    dy = numpy.array([y[1], force(y[0]) / mass])

    return dy


def ralston(t0, y0, t1, f):
    h = t1-t0

    k1 = h*f(t0, y0)
    
    p = y0 + 2.0/3.0 * k1   
    k2 = h*f(t0+2.0/3.0*h, p)

    y1 = y0 + (k1 + 3.0 * k2)/4.0

    return y1



def main():
    t0 = 0.0
    x0 = 2.5e-2
    v0 = 0.0
    
    nsteps = 100000
    h = 0.001

    y0 = numpy.array([x0, v0])
    out = open("sysrk.txt","w")
    
    for i in range(nsteps+1):
        s = str(t0) + ' '+ str(y0[0]) + ' ' + str(y0[1]) + '\n'
        out.write(s)
	
        y1 = ralston(t0, y0, t0+h, f)
	
        t0 += h
        y0 = y1


    out.close()

    
main()
