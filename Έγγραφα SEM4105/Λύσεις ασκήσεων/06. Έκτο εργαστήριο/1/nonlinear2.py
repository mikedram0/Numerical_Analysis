import cmath


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


def g(x,y):
  return 4.0 * x*x - pow(y,3) + 28.0


def h(x,y):
  return 3.0 * pow(x,3) + 4.0 * y*y - 145.0



def f(z):
    x = z.real
    y = z.imag

    re = g(x,y)
    im = h(x,y)

    return complex(re,im)





def main():
    
    z1 = complex(0.0,0.0)
    z2 = complex(1.0,1.0)

    z = secant(z1, z2, 1e-8, f)

    print(z.real, z.imag)


main()
