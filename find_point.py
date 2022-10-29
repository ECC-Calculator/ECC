#montgomery curve : By^2 = x^3 + Ax^2 + x

import gmpy2, random, sys
import graph_points as graph

x_coordinates = []
y_coordinates = []

#gets us the next prime if number isn't prime
def getPrime(number):
    if gmpy2.is_prime(number):
        return number
    else:
        next_prime = gmpy2.next_prime(number)
        print("{} isn't a prime so we consider the next prime {} for calculations".format(number, next_prime))
        return next_prime

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

#y^2 = m (mod p); here m is the function of x
def findM(a, b, x, p):
    return gmpy2.f_mod((x*x*x + a*x*x + x) * (gmpy2.invert(b, p)), p)

def findQRForSW(a, b, x, p):
    return gmpy2.f_mod((x*x*x + a*x + b), p)

def legendre(a, p):
    return modular_pow(a, (p - 1) // 2, p)


def power_modulo(a: int, b: int, n: int) -> int:
    """ Computes a ^ b mod n """
    result = 1
    while b != 0:
        if b % 2 == 1:
            # b odd
            result = (result * a) % n
        a = (a * a) % n
        b >>= 1
    return result


def extended_gcd(a: int, b: int) -> (int, int, int):
    # optional check
    if a == 0:
        return b, 0, 1

    # without this check the first iteration will divide by zero
    if b == 0:
        return a, 1, 0

    un_prev = 1
    vn_prev = 0
    un_cur = 0
    vn_cur = 1

    while True:
        qn = a // b
        new_r = a % b
        a = b
        b = new_r

        if b == 0:
            return a, un_cur, vn_cur

        # Update coefficients
        un_new = un_prev - qn * un_cur
        vn_new = vn_prev - qn * vn_cur

        # Shift coefficients
        un_prev = un_cur
        vn_prev = vn_cur
        un_cur = un_new
        vn_cur = vn_new


def inverse_modulo(a: int, n: int) -> int:
    _, b, _ = extended_gcd(a, n)
    return b % n


def legendre_symbol(a: int, p: int, /) -> int:
    return power_modulo(a, (p - 1) >> 1, p)


def _choose_b(p: int, /) -> int:
    b = 2
    while legendre_symbol(b, p) == 1:
        b = random.randrange(2, p)
    return b


def _tonelli_shanks_recursive(a: int, k: int, p: int, b: int, b_inverse: int, /):
    """
    Computes a square root of a modulo prime p
    :param a: the number to take the square root of
    :param k: positive integer, such that a^m = 1 (mod p) where m = (p-1)/(2^k)
    :param p: odd prime p modulo which we are working
    :param b: an arbitrary non-square modulo p
    :param b_inverse: the inverse of b modulo p, i.e., b * b_inverse = 1 (mod p)
    :return: one of the square roots of a modulo p (the other can be obtained via negation modulo p)
    """

    m = (p - 1) >> k
    a_m = 1

    while m % 2 == 0 and a_m == 1:
        m >>= 1
        k += 1
        a_m = power_modulo(a, m, p)

    if a_m == p - 1:
        # a^m = -1 (mod p)
        b_power = 1 << (k - 1)
        b_power_half = 1 << (k - 2)
        a_next = (a * power_modulo(b, b_power, p)) % p
        a_next_root = _tonelli_shanks_recursive(a_next, k, p, b, b_inverse)
        a_root = a_next_root * power_modulo(b_inverse, b_power_half, p)
        return a_root % p

    # we now handle the case when m is odd
    # this case is easy, a^((m+1)/2) is a square root of a
    return power_modulo(a, (m + 1) >> 1, p)


def tonelli_shanks(a: int, p: int, /):
    """
    Computes a square root of a modulo prime p
    :param a: the number to take the square root of
    :param p: odd prime p modulo which we are working
    :return: one of the square roots of a modulo p (the other can be obtained via negation modulo p)
    """

    if legendre_symbol(a, p) != 1:
      # a is not not a square modulo p
      return None

    b = _choose_b(p)
    b_inverse = inverse_modulo(b, p)
    return _tonelli_shanks_recursive(a, 1, p, b, b_inverse)


def addPoint(x, m, p, A, B):
    if gmpy2.f_mod(p, 4) == 3:			
      y = euler(m, p)
    else:			
      y = tonelli_shanks(m, p)

    print("({}, {})  ({}, {})".format(x, y, x, p-y))
    x_coordinates.extend([int(x), int(x)])
    y_coordinates.extend([int(y), int(p-y)])

def findPoint(a, b, p):
    totalPoints = 0
    p = getPrime(p)
    print("Prime field is:", p)
    print("")
    for i in range(p):
        m = findM(a, b, i, p)
        
        quadraticResidue = legendre(m, p)
        if quadraticResidue == 1:
            addPoint(i, m, p, a, b)
            totalPoints += 2
        elif quadraticResidue == 0:
            print("({}, {})".format(i, 0))
            x_coordinates.append(i)
            y_coordinates.append(0)
            totalPoints += 1

    print("\nNumber of total points (including 0 points) is:", totalPoints)
    print("---------------------------------------------")
    graph.plot(x_coordinates, y_coordinates)

n = len(sys.argv)
if n != 4 or gmpy2.mpz(sys.argv[1]) == 2 or gmpy2.mpz(sys.argv[1]) == -2 or gmpy2.mpz(sys.argv[2]) == 0:
    print("INVALID INPUT! Please try again")
    print("Unable to process without 4 arguments. Provide input as per below format.")
    print("format: python <file_name>.py <A> <B> <p> where A != + 2 and A != -2 and B != 0")    
else:
    print("---------------------------------------------")
    print("Montgomery Curve: By^2 = x^3 + Ax^2 + x\n")
    print("A =", sys.argv[1])
    print("B =", sys.argv[2])
    print("User given prime number =", sys.argv[3])
    print("---------------------------------------------")
    findPoint(gmpy2.mpz(sys.argv[1]), gmpy2.mpz(sys.argv[2]), gmpy2.mpz(sys.argv[3]))