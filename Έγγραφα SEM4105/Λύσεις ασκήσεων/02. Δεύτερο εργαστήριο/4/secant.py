import math

def secant(x1,x2, tol, f):
    f1 = f(x1)
    f2 = f(x2)

    while abs(f2) > tol :
        x = x2 - f2 * (x2 - x1) / (f2 - f1)

        x1 = x2
        f1 = f2

        x2 = x
        f2 = f(x)

    return x


def f(x):
    return 3.0 * math.log(x) + 5.0


def main():
    tol = 1e-8

#     Initial approximations
    a = 0.1
    b = 0.2
    
    x = secant(a,b,tol, f)

    print ("A root is approximately ", x)
    print ("The function value is ", f(x))


main()
