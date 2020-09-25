from convolution import conv_mul_mod, conv_inv, add_mod
from modulo import mod
import numpy as np


def T(d1, d2, n):
    ones = [1] * d1
    nones = [-1] * d2
    zeroes = [0] * (n - d1 - d2)
    arr = np.array(ones + zeroes + nones)
    np.random.shuffle(arr)
    return arr


def center_lift(a, q):
    a = mod(a, q)
    a[a > q / 2] = a[a > q / 2] - q
    return np.array(a)


def public_key(f, g, n, q):
    return conv_mul_mod(conv_inv(f, n, q), g, n, q)


def init(n, p, q, d):
    g = T(d, d, n)
    f = T(d + 1, d, n)
    Fp = conv_inv(f, n, p)
    Fq = conv_inv(f, n, q)
    while Fp == -1 or Fq == -1:
        f = T(d + 1, d, n)
        Fp = conv_inv(f, n, p)
        Fq = conv_inv(f, n, q)
    h = public_key(f, g, n, q)
    return f, g, h


def encrypt(m, r, h, n, p, q):
    mul = conv_mul_mod(h, r, n, q) * p
    return add_mod(mul, m, q)


def decrypt(e, f, n, p, q, d):
    a = conv_mul_mod(f, e, n, q)
    a = center_lift(a, q)
    Fp = conv_inv(f, n, p)
    b = conv_mul_mod(Fp, a, n, p)
    return center_lift(b, p)


def example():
    n = 7
    p = 3
    q = 41
    d = 2
    f = np.array([-1, 0, 1, 1, -1, 0, 1])
    g = np.array([0, -1, -1, 0, 1, 0, 1])
    h = public_key(f, g, n, q)
    print(h)
    m = np.array([1, -1, 1, 1, 0, -1])
    r = np.array([-1, 1, 0, 0, 0, -1, 1])
    e = encrypt(m, r, h, n, p, q)
    print(e)
    m = decrypt(e, f, n, p, q, d)
    print(m)


if __name__ == '__main__':
    # print(T(2, 2, 5))
    # q = 7
    # a = np.array([5, 3, -6, 2, 4])
    # print(center_lift(a, q))
    example()
