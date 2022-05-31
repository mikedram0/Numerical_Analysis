import math

def difeq(x,y,z):
    dy=list((0,0,0,0))
    dz=list((0,0,0,0))


    dy[0] = y + z*z - x**3
    dz[0] = y**3 + z + math.cos(x)

    dy[1] =-3.0 * x*x + dy[0] + 2.0 * z * dz[0]			 
    dz[1] = 3.0 * y*y * dy[0] + dz[0] - math.sin(x)
 
    dy[2] = -6.0 * x + dy[1] + 2.0 * (dz[0]*dz[0] + z * dz[1])
    dz[2] = 3.0 * y * (2.0 * dy[0]*dy[0] + y * dy[1]) + dz[1] - math.cos(x)

    dy[3] = -6.0 + dy[2] + 2.0 * (3.0 * dz[0] * dz[1] + z * dz[2])
    dz[3] = dz[2] + math.sin(x) + 6.0 * dy[0] * dy[0] * dy[0] + 3.0 * y * (6.0 * dy[0] * dy[1] + y * dy[2])

    return dy,dz


def taylorstep(x0,y0,z0,x1):
    dy,dz = difeq(x0,y0,z0)
    y1 = y0
    z1 = z0
    h = x1-x0
    term = h

    for i in range(len(dy)):
        y1 += term * dy[i]
        z1 += term * dz[i]
        term *= h / (i+1)

    return y1,z1
    


def taylor(a,b,h,ya,za):
    xold = a
    yold = ya
    zold = za

    while xold < b:
        xnew = xold + h
        if xnew > b:
            xnew = b
	
        ynew,znew = taylorstep(xold, yold, zold, xnew)

        xold = xnew
        yold = ynew
        zold = znew

        print(xnew,ynew,znew)
        
    return ynew,znew



def main():
    a = 0.0
    b = 1.0
    h = 0.01
    
    yinit = 0.3
    zinit = 0.1

    taylor(a,b,h,yinit,zinit)



main()
