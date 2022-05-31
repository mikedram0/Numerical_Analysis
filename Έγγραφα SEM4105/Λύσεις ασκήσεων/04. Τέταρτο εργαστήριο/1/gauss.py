import numpy


def triang(a,b):
    "Τριγωνοποίηση"

    n = len(b)

    for k in range(n-1):
        for i in range(k+1,n):
            ell = -a[i,k]/a[k,k]

            for j in range(k,n):
                a[i,j] += ell * a[k,j]

            b[i] += ell * b[k]
        


def backsub(a,b,x):
    "Οπισθοδρόμηση"
    
    n = len(b)
    for i in range(n-1,-1,-1):
        x[i] = b[i]
        for j in range(i+1,n):
            x[i] -= a[i,j] * x[j]
        
        x[i] /= a[i,i]


        

def main():
    f = open("gauss1.txt", "r")
    n = int(f.readline())

    x = numpy.ndarray(shape=(n), dtype=numpy.float64)

    a = numpy.fromfile(f, dtype=numpy.float64, count=n*n,sep=" ")
    a = numpy.reshape(a,newshape=(n,n))
    
    b = numpy.fromfile(f, dtype=numpy.float64, count=n,  sep=" ")

    f.close()
    
    triang(a,b)
    backsub(a,b,x)
    
    print (x)

main()
