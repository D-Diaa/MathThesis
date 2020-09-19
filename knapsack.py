from modulo import *


def knapsack_super_inc(ms, s):
    ret = [0] * len(ms)
    for i in range(len(ms) - 1, -1, -1):
        if s >= ms[i]:
            ret[i] = 1
            s -= ms[i]
    return ret


def check_solve_knapsack_superinc(ms, s):
    ans = knapsack_super_inc(ms, s)
    sum = 0
    for i in range(len(ans)):
        sum += ans[i] * ms[i]
    if s == sum:
        print(f"{ans} - Answer is correct")
    else:
        is_super = True
        for i in range(len(ms) - 1):
            if 2 * ms[i] > ms[i + 1]:
                is_super = False
        if not is_super:
            print(f"{ans} - Incorrect because sequence is not superincreasing")
        else:
            print(f"{ans} - Sum is not a solution")


def get_private_sequence(public_sequence, multiplier, modulus):
    _inv = inverse_mod(multiplier, modulus)
    r = []
    for pub in public_sequence:
        r.append(mod(_inv * pub, modulus))
    return r


def Q2():
    check_solve_knapsack_superinc([3, 7, 19, 43, 89, 195], 260)
    check_solve_knapsack_superinc([5, 11, 25, 61, 125, 261], 408)
    check_solve_knapsack_superinc([2, 5, 12, 28, 60, 131, 257], 334)
    check_solve_knapsack_superinc([4, 12, 15, 36, 75, 162], 214)


def Q3():
    A = 4392
    B = 8387
    S = 26560
    S_bar = mod(inverse_mod(A, B)*S, B)
    seq = get_private_sequence([5186, 2779, 5955, 2307, 6599, 6771, 6296, 7306, 4115, 7039], A, B)
    check_solve_knapsack_superinc(seq, S_bar)


if __name__ == '__main__':
    Q2()
    Q3()
