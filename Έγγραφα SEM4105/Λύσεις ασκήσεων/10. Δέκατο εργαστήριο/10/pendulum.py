import numpy
import math


"""
  The system is 
  D theta / D t  = z
  D z / D t      = -sin(theta)

  theta -> y[0], z-> y[1].  D theta / D t -> f[0], D z / D t -> f[1]. t -> x
"""

def f(x,y):
    return numpy.array([y[1], -math.sin(y[0])])
                



def rk4(x,y,h,f):
                    
    k1 = h * f(x,y)
                    
    p = y + k1 / 2.0
    k2 = h * f(x+h/2.0,p)
                        
    p = y + k2 / 2.0
    k3 = h * f(x+h/2.0,p)

    p = y + k3
    k4 = h * f(x+h,p)

    ynew = y + (k1 + 2.0 * (k2 + k3) + k4) / 6.0
    xnew = x + h

    return xnew,ynew



"""
  solution of theta'' = - theta:   theta = A cos(t) + B sin(t)
  theta(0) = 45 deg. , theta'(0) = 0 => 
"""
def approx(t):
    a = 45.0/180.0 * math.pi
    return a * numpy.array([math.cos(t), -math.sin(t)]);
                    

def main():
    tinit = 0.0
    tfin = 10.0
    niter = 100
    h = (tfin-tinit)/niter

    yinit = numpy.array([45.0 / 180.0 * math.pi, 0.0])

    t = tinit
    y = yinit
   
    for k in range(niter):
        t,y = rk4(t,y,h,f)

        yapprox = approx(t)

        print(t, y, yapprox)


main()
