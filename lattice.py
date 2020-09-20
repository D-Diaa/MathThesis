import numpy as np
import scipy.linalg as lin
from utils import first_divisor
from vector_spaces import gram_schmidt
from math import pi, e, sqrt
from random import randint

EPS = 1e-9


def eps_eq(u, v):
    return u - EPS <= v <= u + EPS


def is_integer_entries(mat):
    ints = np.rint(mat)
    return np.allclose(ints, mat)


def in_GL(mat):
    det = lin.det(mat)
    return (eps_eq(det, 1) or eps_eq(det, -1)) and is_integer_entries(mat)


def GL_inverse(mat):
    if in_GL(mat):
        print(lin.inv(mat))
    else:
        print("Not in GL")


def gaussian_heuristic(dim, mat=None, det=None, nthroot_det=None):
    sigma = sqrt(dim / (2 * pi * e))
    if nthroot_det is None:
        if det is None:
            if mat is None:
                return -1
            else:
                det = lin.det(mat)
        nthroot_det = det ** (1 / dim)
    return sigma * nthroot_det


def gen_good_basis(dim, min_val=-100, max_val=100, min_h=0.8):
    h = 0
    mat = None
    while h < min_h:
        mat = np.array(np.random.randint(min_val, max_val, (dim, dim)), dtype='float64')
        h = hadamard_ratio(mat)
        if h < min_h:
            mat = gram_schmidt(mat)
            mat = np.rint(mat)
            h = hadamard_ratio(mat)
        if h != h:
            h = 0
    # print(f"found basis with h = {h}")
    return mat


def randomize_rows(mat):
    dim = len(mat)
    for i in range(dim):
        swapInd = randint(i, dim - 1)
        mat[[i, swapInd]] = mat[[swapInd, i]]
    return mat


def gen_rand_GLN_single(dim, max_val=100):
    max = int(sqrt(2 * max_val / dim))
    L = np.identity(dim)
    U = np.identity(dim)
    for i in range(0, dim):
        for j in range(0, i):
            L[i][j] = randint(1, max)
            U[j][i] = randint(1, max)
    ret = np.matmul(L, U)
    ret = randomize_rows(ret)
    ret = np.transpose(ret)
    ret = randomize_rows(ret)
    return ret


def gen_rand_GL2(samples, max_val=100):
    res = []
    while len(res) <= samples:
        b = randint(1, max_val)
        c = randint(1, max_val)
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


def babai_cvp(bases, w):
    bases = np.transpose(bases)
    _inv = lin.inv(bases)
    t = np.matmul(_inv, w)
    a = np.rint(t)
    return np.matmul(bases, a)


def hadamard_ratio(bases):
    dim = len(bases)
    det = abs(lin.det(bases))
    denom = 1
    for base in bases:
        denom *= lin.norm(base)
    frac = det / denom
    return frac ** (1 / dim)


def babai_svp(bases):
    pass


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


def Q16():
    print(gaussian_heuristic(251, nthroot_det=2 ** (2251.58 / 251)))


def Q17():
    basis = [
        [213, -437],  # v1
        [312, 105]  # v2
    ]
    w = [43127, 11349]

    v = babai_cvp(basis, w)
    h = hadamard_ratio(basis)
    d = lin.norm(v - w)
    print(v)
    print(d)
    print(h)
    ans = "Very good basis" if h >= 0.9 else "basis not very orthogonal"
    print(ans)

    new_basis = [
        [2937, -1555],
        [11223, -5888]
    ]
    M = change_of_basis(basis, new_basis)
    print(M)
    ans = "new basis is also a basis for L" if in_GL(M) else "not a new basis"
    print(ans)

    v = babai_cvp(new_basis, w)
    h = hadamard_ratio(new_basis)
    d = lin.norm(v - w)
    print(v)
    print(d)
    print(h)
    ans = "Very good basis" if h >= 0.9 else "basis not very orthogonal"
    print(ans)


if __name__ == '__main__':
    # Q7()
    # Q11()
    # Q12()
    # Q13()
    # Q16()
    Q17()
