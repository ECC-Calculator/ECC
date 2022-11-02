# Montgomery Curve Equation : by^2 = x^3 + ax^2 + x (mod p)

import gmpy2, sys

# Adding 2 points (x1, y1) & (x2, y2)
def pointAddition(x1, y1, x2, y2, a, b, p):
    if x1 == x2 and y1 == y2 :
        pointDoubling(x1, y1, a, b, p)
    elif x1 == x2 :
        return 0, 0
    else :
        try :
            k = gmpy2.mul(gmpy2.sub(y2, y1) % p, gmpy2.invert(gmpy2.sub(x2, x1) % p, p)) % p
            x3 = gmpy2.sub(gmpy2.sub(gmpy2.sub(gmpy2.mul(b, gmpy2.powmod(k, 2, p)) % p, a) % p, x1) % p, x2) % p
            y3 = gmpy2.sub(gmpy2.sub(gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(2, x1) % p, x2) % p, a) % p, k) % p, gmpy2.mul(b, gmpy2.powmod(k, 3, p)) % p), y1) % p
            
            if x3 < 0:
                x3 = x3 + p
            if y3 < 0:
                y3 = y3 + p
            
            return x3, y3
        except Exception as e: 
            print(e)

# Point Doubling of the point (x, y)
def pointDoubling(x, y, a, b, p):
    if y == 0 :
        return 0, 0
    try :
        k = gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(gmpy2.powmod(x, 2, p), 3) % p, gmpy2.mul(gmpy2.mul(2, a) % p, x) % p) % p, 1) % p, gmpy2.invert(gmpy2.mul(gmpy2.mul(2, b) % p, y) % p, p))
        x3 = gmpy2.sub(gmpy2.sub(gmpy2.sub(gmpy2.mul(b, gmpy2.powmod(k, 2, p)) % p, a) % p, x) % p, x) % p
        y3 = gmpy2.sub(gmpy2.sub(gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(2, x) % p, x) % p, a) % p, k) % p, gmpy2.mul(b, gmpy2.powmod(k, 3, p)) % p), y) % p

        if x3 < 0:
            x3 = x3 + p
        if y3 < 0:
            y3 = y3 + p
        
        return x3, y3
    except Exception as e: 
        print(e)

# Point Multiplication of the point (x, y) by a scalar k
def pointMultiplication(x, y, k, a, b, p):
    x0 = 0
    y0 = 0
    x1 = x
    y1 = y
    idx = gmpy2.bit_length(k)
    while idx >= 0:
        if gmpy2.bit_test(k, idx):
            x0, y0 = pointAddition(x0, y0, x1, y1, a, b, p)
            x1, y1 = pointDoubling(x1, y1, a, b, p)
        else:
            x1, y1 = pointAddition(x0, y0, x1, y1, a, b, p)
            x0, y0 = pointDoubling(x0, y0, a, b, p)
        idx = idx - 1
    
    print("({}, {}) multiplied by {} is = ({}, {})".format(x, y, k, x0, y0))

n = len(sys.argv)
if n != 7 or gmpy2.mpz(sys.argv[1]) == 2 or gmpy2.mpz(sys.argv[1]) == -2 or gmpy2.mpz(sys.argv[2]) == 0:
    print("INVALID INPUT!!! Please try again!")
    print("Unable to process without 7 arguments. Provide input as per below format.")
    print("format: python <file_name>.py <A> <B> <p> <x> <y> <k> where A != 2, A != -2 and B != 0") 
else:
    print("------------------------------------------------------------")
    print("Montgomery Curve: By^2 = x^3 + Ax^2 + x\n")
    print("A =", sys.argv[1])
    print("B =", sys.argv[2])
    print("User given prime number =", sys.argv[3])
    print("------------------------------------------------------------")
    pointMultiplication(gmpy2.mpz(sys.argv[4]), gmpy2.mpz(sys.argv[5]), gmpy2.mpz(sys.argv[6]), gmpy2.mpz(sys.argv[1]), gmpy2.mpz(sys.argv[2]), gmpy2.mpz(sys.argv[3]))