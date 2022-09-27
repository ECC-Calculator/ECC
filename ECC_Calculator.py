import gmpy2

#gets us the next prime if number isn't prime
def getPrime(number):
	if gmpy2.is_prime(number):
		return number
	else:
		return gmpy2.next_prime(number)

#fast modular exponentiation
def modular_pow(base, exponent, modulus):
	result = 1
	while exponent > 0:
		if exponent & 1:
			result = gmpy2.f_mod(result * base, modulus)
		exponent = exponent >> 1
		base = gmpy2.f_mod(base * base, modulus)
	return result

#euler formula to calculate y
def euler(quadraticResidue, prime):
	return modular_pow(quadraticResidue, gmpy2.f_div(prime+1, 4), prime)

def findQR(a, b, x, p):
    return gmpy2.f_mod((x*x*x + a*x*x + x) // b, p)

def findQRForSW(a, b, x, p):
    return gmpy2.f_mod((x*x*x + a*x + b), p)

def legendre(a, p):
    return pow(a, (p - 1) // 2, p)

def tonelli(n, p):
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r

def addPoint(x, m, p):
    if gmpy2.f_mod(p, 4) == 3:
      y = euler(m, p)
    else:
      y = tonelli(m, p)

    print(x, y)
    print(x, p - y)

def findPoint(a, b, p):
    for i in range(p):
        # m = findQR(a, b, i, p)
        m = findQRForSW(a, b, i, p)
        # qr = (pow(m, (p - 1) // 2)) % p
        quadraticResidue = modular_pow(m, gmpy2.f_div(p-1, 2), p)
        if quadraticResidue == 1:
            addPoint(i, m, p)
        elif quadraticResidue == 0:
            print('Single point solution: ', i, 0)

def driver():
    findPoint(2, 3, 13)

driver()