import math

def f(x,y):
    return math.cos(x) - math.sin(y) + x*x


def forward_euler(x0, y0, x1, f):
    return y0 + (x1-x0) * f(x0, y0)


# y1 = y0 + (x1-x0) f(x1, y1)
def g(y1,x0,y0,x1,f):
    return y1 - y0 - (x1-x0) * f(x1,y1)



def backward_euler(x0, y0, x1, f):
    #secant for g(y, x0,y0,x1,f) = 0
    #need two points: one is y1 from forward_euler
    y1 = forward_euler(x0, y0, x1, f)
    #the other is y2, close to y1
    y2 = y1 * 1.1
    
    # y will be the approximation to the unknown (y1)
    while True:
        g1 = g(y1, x0, y0, x1, f)
        g2 = g(y2, x0, y0, x1, f)
	
        y = y2 -  g2 * (y2-y1) / (g2 - g1)
	
        if abs(g(y, x0, y0, x1, f)) < 1e-8:
            break
	
        y1 = y2
        y2 = y
        
    return y







def main():
    a = -1.0
    yinit = 3.0
    b = 1.0
    h = 0.01
    
    x0 = a
    y0f = yinit
    y0b = yinit

    n = round((b-a)/h)

    print ("#X\tforward\t\tbackward\n")
    print(x0, y0f, y0b)
    
    for i in range(n):
	#change h if x0 is beyond b:
        if x0+h>b:
            h = b-x0

        y0f = forward_euler(x0, y0f, x0+h, f)
        y0b = backward_euler(x0, y0b, x0+h, f)
        x0 += h
	
        print(x0, y0f, y0b)
    
    


main()
