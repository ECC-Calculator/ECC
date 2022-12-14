# Montgomery Curve Equation : by^2 = x^3 + ax^2 + x (mod p)

import gmpy2, sys

# Adding 2 points (x1, y1) & (x2, y2)
def pointAddition(x1, y1, x2, y2, a, b, p):
    if x1 == x2 and y1 == y2 :
        print('Both the points are same! Perform point doubling operation instead addition.')
    elif x1 == x2 :
        print("Addition of two points ({}, {}) and ({}, {}) is = ({}, {})".format(x1, y1, x2, y2, 0, 0))
    else :
        try :
            k = gmpy2.mul(gmpy2.sub(y2, y1) % p, gmpy2.invert(gmpy2.sub(x2, x1) % p, p)) % p
            x3 = gmpy2.sub(gmpy2.sub(gmpy2.sub(gmpy2.mul(b, gmpy2.powmod(k, 2, p)) % p, a) % p, x1) % p, x2) % p
            y3 = gmpy2.sub(gmpy2.sub(gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(2, x1) % p, x2) % p, a) % p, k) % p, gmpy2.mul(b, gmpy2.powmod(k, 3, p)) % p), y1) % p
            
            if x3 < 0:
                x3 = x3 + p
            if y3 < 0:
                y3 = y3 + p
            
            print("Addition of two points ({}, {}) and ({}, {}) is = ({}, {})".format(x1, y1, x2, y2, x3, y3))
        except Exception as e: 
            print(e)

n = len(sys.argv)
if n != 8 or gmpy2.mpz(sys.argv[1]) == 2 or gmpy2.mpz(sys.argv[1]) == -2 or gmpy2.mpz(sys.argv[2]) == 0:
    print("INVALID INPUT!!! Please try again!")
    print("Unable to process without 8 arguments. Provide input as per below format.")
    print("format: python <file_name>.py <A> <B> <p> <x1> <y1> <x2> <y2> where A != 2, A != -2 and B != 0") 
else:
    print("------------------------------------------------------------")
    print("Montgomery Curve: By^2 = x^3 + Ax^2 + x\n")
    print("A =", sys.argv[1])
    print("B =", sys.argv[2])
    print("User given prime number =", sys.argv[3])
    print("------------------------------------------------------------")
    pointAddition(gmpy2.mpz(sys.argv[4]), gmpy2.mpz(sys.argv[5]), gmpy2.mpz(sys.argv[6]), gmpy2.mpz(sys.argv[7]), gmpy2.mpz(sys.argv[1]), gmpy2.mpz(sys.argv[2]), gmpy2.mpz(sys.argv[3]))