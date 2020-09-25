from lattice import *


def init(dim, d, good_h=0.8, bad_h=0.1):
    V = gen_good_basis(dim, -d, d, good_h)
    print(f"good h = {hadamard_ratio(V)}")
    h = 1
    while h > bad_h:
        U = gen_rand_GLN_single(dim, d)
        W = np.matmul(U, V)
        h = hadamard_ratio(W)
    print(f"bad h = {h}")
    return U, V, W


def encrypt(W, m, r):
    return np.matmul(m, W) + r


def decrypt(V, W, e):
    v = babai_cvp(V, e)
    _inv = lin.inv(W)
    return np.rint(np.matmul(v, _inv))


def example():
    V = [
        [-97, 19, 19],
        [-36, 30, 86],
        [-184, -64, 78]
    ]
    print(hadamard_ratio(V))
    U = [
        [4327, -15447, 23454],
        [3297, -11770, 17871],
        [5464, -19506, 29617]
    ]
    W = np.matmul(U, V)
    print(hadamard_ratio(W))
    m = [86, -35, -32]
    r = [-4, -3, 2]
    e = encrypt(W, m, r)
    print(e)
    good_m = decrypt(V, W, e)
    print(good_m)
    bad_m = decrypt(W, W, e)
    print(bad_m)


def ggh_experiment(d=300, delta=15, dim=10, good_h=0.99, bad_h=0.1):
    U, V, W = init(dim, d, good_h, bad_h)
    m = np.random.randint(-d, d, (1, dim))[0]
    r = np.random.randint(-delta, delta, (1, dim))[0]
    print(f"original m: {m}")
    enc = encrypt(W, m, r)
    good_m = np.array(decrypt(V, W, enc), dtype='int64')
    print(f"dec_good m: {good_m}")
    bad_m = np.array(decrypt(W, W, enc), dtype='int64')
    print(f"dec_bad m: {bad_m}")


def Q18():
    V = [
        [4, 13],
        [-57, -45]
    ]
    W = [
        [25453, 9091],
        [-16096, -5749]
    ]
    print(lin.det(V))
    print(f"private: {hadamard_ratio(V)}")
    print(f"public: {hadamard_ratio(W)}")

    e = [155340, 55483]
    m = decrypt(V, W, e)
    r = np.array(e) - np.matmul(m, W)
    print(f"m={m}")
    print(f"r={r}")
    print(f"using bad basis={ decrypt(W, W, e)}")


if __name__ == '__main__':
    # example()
    Q18()
