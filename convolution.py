from modulo import mod, inverse_mod
import numpy as np


def deg(a):
    for i in range(len(a) - 1, 0, -1):
        if a[i] != 0:
            return i
    return 0


def reduce(c):
    return c[:deg(c) + 1]


def extend(a, n):
    ac = np.zeros(n)
    m = deg(a)
    ac[:m + 1] = a[:m + 1]
    return ac


def subtract(a, b):
    n = deg(a)
    m = deg(b)
    l = max(n, m)
    ac = np.zeros(l + 1)
    bc = np.zeros(l + 1)
    ac[:n + 1] = a[:n + 1]
    bc[:m + 1] = b[:m + 1]
    c = ac - bc
    return reduce(c)


def subtract_mod(a, b, q):
    return mod(subtract(a, b), q)


def add(a, b):
    n = deg(a)
    m = deg(b)
    l = max(n, m)
    ac = np.zeros(l + 1)
    bc = np.zeros(l + 1)
    ac[:n + 1] = a[:n + 1]
    bc[:m + 1] = b[:m + 1]
    c = ac + bc
    return reduce(c)


def add_mod(a, b, q):
    return mod(add(a, b), q)


def poly_mul(a, b):
    n = deg(a)
    m = deg(b)
    c = np.zeros(n + m + 1)
    ac = np.zeros(n + m + 1)
    bc = np.zeros(n + m + 1)
    ac[:n + 1] = a[:n + 1]
    bc[:m + 1] = b[:m + 1]
    for i in range(deg(a) + 1):
        c += ac[i] * bc
        bc = np.roll(bc, 1)
    return c


def poly_mul_mod(a, b, q):
    n = deg(a)
    m = deg(b)
    c = np.zeros(n + m + 1)
    ac = np.zeros(n + m + 1)
    bc = np.zeros(n + m + 1)
    ac[:n + 1] = a[:n + 1]
    bc[:m + 1] = b[:m + 1]
    for i in range(deg(a) + 1):
        c += ac[i] * bc
        c = mod(c, q)
        bc = np.roll(bc, 1)
    return c


def poly_divmod(a, b, q):
    n = deg(a)
    d = deg(b)
    a = mod(a, q)
    b = mod(b, q)
    b = b[:d + 1]

    r = a[:n + 1]
    e = deg(r)
    k = np.zeros(n + 1)
    while e >= d and not (e == 0 and r[0] == 0):
        temp = np.zeros(e - d + 1)
        temp[e - d] = mod(r[e] * inverse_mod(b[d], q), q)
        k = add_mod(k, temp, q)
        r = subtract_mod(r, poly_mul_mod(temp, b, q), q)
        e = deg(r)

    return k, r


def poly_update_egcd(r0, r1, k, q):
    tmp = subtract_mod(r0, poly_mul_mod(k, r1, q), q)
    r0 = r1
    r1 = tmp
    return r0, r1


def poly_egcd(a, b, q):
    a = mod(a, q)
    b = mod(b, q)
    r0 = a
    r1 = b
    x1 = y0 = np.array([0])
    y1 = x0 = np.array([1])
    while r1.any():
        k, _ = poly_divmod(r0, r1, q)
        r0, r1 = poly_update_egcd(r0, r1, k, q)
        x0, x1 = poly_update_egcd(x0, x1, k, q)
        y0, y1 = poly_update_egcd(y0, y1, k, q)
    return reduce(r0), reduce(x0), reduce(y0)


def conv_mul(a, b, n):
    if len(a) < n:
        a = extend(a, n)
    if len(b) < n:
        b = extend(b, n)
    c = np.zeros(n)
    for i in range(n):
        for j in range(n):
            k = mod(i + j, n)
            c[k] += a[i] * b[k - i]
    return c


def conv_mul_mod(a, b, n, q):
    if len(a) < n:
        a = extend(a, n)
    if len(b) < n:
        b = extend(b, n)
    c = np.zeros(n)
    for i in range(n):
        for j in range(n):
            k = mod(i + j, n)
            c[k] = mod(c[k] + a[i] * b[k - i], q)
    return c


def conv_inv(a, n, q):
    b = np.zeros(n + 1)
    b[0] = -1
    b[n] = 1
    g, x, y = poly_egcd(b, a, q)
    if deg(g) > 0 or g[0] != 1:
        return -1
    return y


def example():
    n = 5
    a = [1, -2, 0, 4, -1]
    b = [3, 4, -2, 5, 2]
    print(conv_mul(a, b, n))


def test_poly():
    a = np.array([0, 1, 0, 2, 0, 3, 0])
    b = np.array([1, 0, 2, 0, 3, 0, 4])
    print(deg(a))
    print(deg(b))
    print(poly_mul(a, b))
    print(poly_mul_mod(a, b, 10))
    q = 13
    a = np.array([-1, 0, 0, 0, 0, 1])
    b = np.array([-3, 2, 0, 1])
    g, x, y = poly_egcd(a, b, q)
    print(g)
    print(x)
    print(y)
    n = 5
    q = 2
    a = np.array([1, 1, 0, 0, 1])
    _inv = conv_inv(a, n, q)
    print(_inv)
    print(f"{a}*{_inv}={conv_mul_mod(a, _inv, n, q)}")


def Q23():
    n = 3
    q = 7
    a = [1, 1, 0]
    b = [-5, 4, 2]
    print(conv_mul_mod(a, b, n, q))
    n = 10
    q = 2
    a = [0] * n
    a[2] = a[5] = a[7] = a[8] = a[9] = 1
    b = [1] * n
    b[2] = b[6] = 0
    print(conv_mul_mod(a, b, n, q))


def Q25():
    n = 5
    q = 3
    b = np.array([1, 0, 1, -1])
    print(conv_inv(b, n, q))


if __name__ == "__main__":
    # example()
    # Q23()
    # test_poly()
    Q25()
