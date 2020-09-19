from modulo import *


def public_key(f, g, q):
    return mod(inverse_mod(f, q) * g, q)


def encrypt(m, r, h, q):
    return mod(r * h + m, q)


def decrypt(e, f, g, q):
    a = mod(f * e, q)
    f_inv_g = inverse_mod(f, g)
    m = mod(f_inv_g * a, g)
    return m


def Q1():
    f = 19928
    g = 18643
    q = 918293817
    h = public_key(f, g, q)
    m = decrypt(619168806, f, g, q)
    e = encrypt(10220, 19564, h, q)
    print(f"a - h = {h}")
    print(f"b - m = {m}")
    print(f"c - e = {e}")


if __name__ == '__main__':
    Q1()
