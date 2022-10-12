import gmpy2

# Montgomery Curve Equation : by^2 = x^3 + ax^2 + x (mod p)

# Adding 2 points (x1, y1) & (x2, y2)
def pointAddition(x1, y1, x2, y2, a, b, p):
    k = gmpy2.mul(gmpy2.sub(y2, y1) % p, gmpy2.invert(gmpy2.sub(x2, x1) % p, p)) % p
    x3 = gmpy2.sub(gmpy2.sub(gmpy2.sub(gmpy2.mul(b, gmpy2.powmod(k, 2, p)) % p, a) % p, x1) % p, x2) % p
    y3 = gmpy2.sub(gmpy2.sub(gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(2, x1) % p, x2) % p, a) % p, k) % p, gmpy2.mul(b, gmpy2.powmod(k, 3, p)) % p), y1) % p
    
    if x3 < 0:
        x3 = x3 + p
    if y3 < 0:
        y3 = y3 + p
    
    print('Performing Point Addition...')
    print('Resultant Point :', x3, y3)

# Point Doubling of the point (x, y)
def pointDoubling(x, y, a, b, p):
    k = gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(gmpy2.powmod(x, 2, p), 3) % p, gmpy2.mul(gmpy2.mul(2, a) % p, x) % p) % p, 1) % p, gmpy2.invert(gmpy2.mul(gmpy2.mul(2, b) % p, y) % p, p))
    x3 = gmpy2.sub(gmpy2.sub(gmpy2.sub(gmpy2.mul(b, gmpy2.powmod(k, 2, p)) % p, a) % p, x) % p, x) % p
    y3 = gmpy2.sub(gmpy2.sub(gmpy2.mul(gmpy2.add(gmpy2.add(gmpy2.mul(2, x) % p, x) % p, a) % p, k) % p, gmpy2.mul(b, gmpy2.powmod(k, 3, p)) % p), y) % p

    if x3 < 0:
        x3 = x3 + p
    if y3 < 0:
        y3 = y3 + p
    
    print('Performing Point Doubling...')
    print('Resultant Point :', x3, y3)

def driver(a, b, p, x1, y1, x2, y2):
    if x1 == x2 and y1 == y2:
        pointDoubling(x1, y1, a, b, p)
    else:
        pointAddition(x1, y1, x2, y2, a, b, p)

driver(2, 3, 11, 1, 2, 4, 5)