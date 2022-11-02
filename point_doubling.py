# Montgomery Curve Equation : by^2 = x^3 + ax^2 + x (mod p)

import gmpy2, sys

# Point Doubling of the point (x, y)
def pointDoubling(x, y, a, b, p):
    try :
        k = gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(gmpy2.powmod(x, 2, p), 3) % p, gmpy2.mul(gmpy2.mul(2, a) % p, x) % p) % p, 1) % p, gmpy2.invert(gmpy2.mul(gmpy2.mul(2, b) % p, y) % p, p))
        x3 = gmpy2.sub(gmpy2.sub(gmpy2.sub(gmpy2.mul(b, gmpy2.powmod(k, 2, p)) % p, a) % p, x) % p, x) % p
        y3 = gmpy2.sub(gmpy2.sub(gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(2, x) % p, x) % p, a) % p, k) % p, gmpy2.mul(b, gmpy2.powmod(k, 3, p)) % p), y) % p

        if x3 < 0:
            x3 = x3 + p
        if y3 < 0:
            y3 = y3 + p
        
        print('Performing Point Doubling on the point ({}, {})'.format(x, y))
        print('Resultant Point : ({}, {})'.format(x3, y3))
    except Exception as e: 
        print(e)

n = len(sys.argv)
if n != 6 or gmpy2.mpz(sys.argv[1]) == 2 or gmpy2.mpz(sys.argv[1]) == -2 or gmpy2.mpz(sys.argv[2]) == 0:
    print("INVALID INPUT!!! Please try again!")
    print("Unable to process without 8 arguments. Provide input as per below format.")
    print("format: python <file_name>.py <A> <B> <p> <x> <y> where A != 2, A != -2 and B != 0") 
else:
    print("------------------------------------------------------------")
    print("Montgomery Curve: By^2 = x^3 + Ax^2 + x\n")
    print("A =", sys.argv[1])
    print("B =", sys.argv[2])
    print("User given prime number =", sys.argv[3])
    print("------------------------------------------------------------")
    pointDoubling(gmpy2.mpz(sys.argv[4]), gmpy2.mpz(sys.argv[5]), gmpy2.mpz(sys.argv[1]), gmpy2.mpz(sys.argv[2]), gmpy2.mpz(sys.argv[3]))