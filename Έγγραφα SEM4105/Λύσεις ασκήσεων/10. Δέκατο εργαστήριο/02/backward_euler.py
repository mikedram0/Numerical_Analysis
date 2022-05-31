import math


def correct(x):
    return (math.exp(-50.0*x)+ 2500.0*math.cos(x) + 50.0*math.sin(x))/2501.0

def f(x,y):
    return (math.cos(x) - y)/0.02


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
    a = 0.0
    yinit = 1.0
    b = 0.2
    h = 1e-3
    
    x0 = a
    y0 = yinit

    n = round((b-a)/h)

    print(x0, y0, correct(x0))
    
    for i in range(n):
	#change h if x0 is beyond b:
        if x0+h>b:
            h = b-x0
		
        y0 = backward_euler(x0, y0, x0+h, f)
        x0 += h
	
        print(x0, y0, correct(x0))
    
    


main()
