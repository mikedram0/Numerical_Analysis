import math

def difeq(x,y):
    dy=list((0,0,0,0,0))

    dy[0] = math.cos(x) - math.sin(y) + x*x
    dy[1] = 2.0 * x - math.sin(x) - math.cos(y) * dy[0]    
    dy[2] = 2.0 - math.cos(x) -  math.cos(y) * dy[1] + math.sin(y) * dy[0]**2
    
    dy[3] =  math.sin(x) + 3.0 *  math.sin(y) * dy[0] * dy[1] + math.cos(y) * (dy[0]**3 - dy[2])
	
    dy[4] =  math.cos(x) + math.cos(y) * (6.0 * dy[0]**2 * dy[1] - dy[3]) + (3.0 * dy[1]**2 + dy[0] * (4.0 * dy[2] - dy[0]**3)) * math.sin(y)

    return dy


def taylorstep(x0,y0,x1):
    dy = difeq(x0,y0)
    y1 =y0
    h = x1-x0
    term = h

    for i in range(len(dy)):
        y1 += term * dy[i]
        term *= h / (i+2)

    return y1
    



def taylor(a,b,h,ya):
    xold = a
    yold = ya

    while xold < b:
        xnew = xold + h
        if xnew > b:
            xnew = b
	
        ynew = taylorstep(xold, yold, xnew)

        xold = xnew
        yold = ynew
    
    return ynew



def main():
    xa = -1.0
    xb = 1.0
    h = 0.01
    ya = 3.0
    yb = taylor(xa, xb, h, ya)
    
    print(xb, yb)




main()
