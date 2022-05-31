import math

def f(x,y):
    return math.cos(x)-x*math.sin(x)


def forward_euler(x0, y0, x1, f):
    return y0 + (x1-x0) * f(x0, y0)


		     
def correct(x):
    return 2.0 + x*math.cos(x)


def main():
    a = 0.0
    yinit = 2.0
    b = 3.0
    
    h = 0.01
    
    x0 = a
    y0 = yinit

    print(x0, y0, correct(x0))

    while x0 < b:
        if x0+h > b:
            h = b-x0

        y0 = forward_euler(x0, y0, x0+h, f)
        x0 = x0+h
        print(x0, y0, correct(x0))
      


main()
