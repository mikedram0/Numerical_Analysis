import math



def g(x):
    return math.sin(x)


def trapez(f, a,b,n):
    
    h = (b-a)/n
    s = (f(a)+f(b))/2.0

    for i in range(1,n):
        s += f(a+i*h)
    
    return s*h


def main():
    n = 2
    correct = 2.0
    for k in range(1,10):
        integr = trapez(g, 0.0, math.pi, n)
        print(n, integr, correct - integr)
        n *= 2




main()
