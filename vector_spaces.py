import numpy as np
import numpy.linalg as lin


def gram_schmidt(bases):
    for i in range(1, len(bases)):
        for j in range(i):
            mu = np.dot(bases[i], bases[j]) / np.dot(bases[j], bases[j])
            bases[i] -= mu * bases[j]
    return bases


def Q5():
    B = [
        [1, 3, 2],
        [2, -1, 3],
        [1, 0, 2]
    ]
    B_bar = [
        [-1, 0, 2],
        [3, 1, -1],
        [1, 0, 1]
    ]
    M = np.matmul(B, lin.inv(B_bar))
    print(M)


def Q6():
    v1 = np.array([1, 3, 2], dtype='float64')
    v2 = np.array([4, 1, -2], dtype='float64')
    v3 = np.array([-2, 1, 3], dtype='float64')
    print(gram_schmidt([v1, v2, v3]))

    v1 = np.array([4, 1, 3, -1], dtype='float64')
    v2 = np.array([2, 1, -3, 4], dtype='float64')
    v3 = np.array([1, 0, -2, 7], dtype='float64')
    print(gram_schmidt([v1, v2, v3]))


if __name__ == '__main__':
    Q5()
    Q6()
