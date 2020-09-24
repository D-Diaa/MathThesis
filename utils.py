def first_divisor(n):
    i = 2
    while i * i <= n:
        if n % i == 0:
            return i
        i += 1
    return n


def divisors(n):
    small = []
    large = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            small.append(i)
            large.append(n // i)
        i += 1
    large.reverse()
    small.extend(large)
    return small
