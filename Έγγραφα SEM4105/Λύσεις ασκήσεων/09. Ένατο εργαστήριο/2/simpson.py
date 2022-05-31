import math

def g(x):
    return math.sin(x)


def simpson(f, a,b,n):
    if n%2 != 0 : print("error\n")
    
    h = (b-a)/n
    s = f(a)+f(b)

    s1 = 0.0
    for i in range(1,n,2):
        s1 += f(a+i*h)
    
    s += 4.0 * s1

    s2 = 0.0
    for i in range(2,n,2):
        s2 += f(a+i*h)

    s += 2.0 * s2
    
    return s*h/3.0


def main():
    n = 2
    correct = 2.0
    for k in range(1,10):
        integr = simpson(g, 0.0, math.pi, n)
        print(n, integr, correct - integr)
        n *= 2




main()

