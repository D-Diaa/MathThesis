import numpy as np
import numpy.linalg as lin
from utils import first_divisor
from random import randint

EPS = 1e-9


def eps_eq(u, v):
    return u - EPS <= v <= u + EPS


def in_GL(mat):
    det = lin.det(mat)
    return eps_eq(det, 1) or eps_eq(det, -1)


def GL_inverse(mat):
    if in_GL(mat):
        print(lin.inv(mat))
    else:
        print("Not in GL")


def gen_rand_GL2(samples, max=100):
    res = []
    while len(res) <= samples:
        b = randint(1, max)
        c = randint(1, max)
        div = first_divisor(b * c + 1)
        if div != b * c + 1:
            a = div
            d = (b * c + 1) // a
            res.append([
                [a, b],
                [c, d]
            ])
        div = first_divisor(b * c - 1)
        if div != b * c - 1:
            a = div
            d = (b * c + 1) // a
            res.append([
                [a, b],
                [c, d]
            ])
    return res


def check_commutative(mat1, mat2):
    if mat1 == mat2:
        return True
    res1 = np.matmul(mat1, mat2)
    res2 = np.matmul(mat2, mat1)
    return (res1 == res2).all()


def change_of_basis(A, B):
    return np.matmul(B, lin.inv(A))


def Q7():
    mat = np.array([
        [1, 3, -2],
        [2, 1, 0],
        [-1, 2, 5]
    ])
    vol = abs(lin.det(mat))
    print(vol)


def Q11():
    sample = gen_rand_GL2(3, 3)
    checks = [(-1, -1) if check_commutative(a, b) else (a, b) for a in sample for b in sample]
    for check in checks:
        if check != (-1, -1):
            print(check[0])
            print(check[1])
            return


def Q12():
    GL_inverse([
        [3, 1],
        [2, 2]
    ])
    GL_inverse([
        [3, -2],
        [2, -1]
    ])
    GL_inverse([
        [3, 2, 2],
        [2, 1, 2],
        [-1, 3, 1]
    ])
    GL_inverse([
        [-3, -1, 2],
        [1, -3, -1],
        [3, 0, -2]
    ])


def Q13():
    B = [
        [3, 1, -2],
        [1, -3, 5],
        [4, 2, 1]
    ]
    B_1 = [
        [5, 13, -13],
        [0, -4, 2],
        [-7, -13, 18]
    ]
    B_2 = [
        [4, -2, 3],
        [6, 6, -6],
        [-2, -4, 7]
    ]
    M1 = change_of_basis(B, B_1)
    det1 = lin.det(M1)
    M2 = change_of_basis(B, B_2)
    det2 = lin.det(M2)

    if eps_eq(det1, 1) or eps_eq(det1, -1):
        print("B1 is a basis")
    else:
        print("Not a basis")
    if eps_eq(det2, 1) or eps_eq(det2, -1):
        print("B2 is a basis")
    else:
        print("Not a basis")


if __name__ == '__main__':
    Q7()
    Q11()
    Q12()
    Q13()
