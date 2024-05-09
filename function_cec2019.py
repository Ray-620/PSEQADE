import math
import numpy as np

def storns_chebyshev_polynomial_fitting(x):#D=9 search range{-8192,8192}
    u = 0
    v = 0
    wk = 0
    p = 0
    d= 72.661
    m = 32*len(x)
    for i in range(len(x)):
        u += x[i]*((1.2)**(len(x)-i-1))
        v += x[i] * ((-1.2) ** (len(x) - i - 1))
        for k in range(m+1):
            wk += x[i] * (((2*k)/(m)-1) ** (len(x) - i - 1))
            if wk > 1:
                p += (wk-1)**2
            elif wk  < 1:
                p += (wk+1)**2
            else:
                p += 0
    if u < d:
        p1 = (u-d)**2
    else:
        p1 = 0
    if v < d:
        p2 = (v-d)**2
    else:
        p2 = 0
    return p1 + p2 +p

def ackley(x):
    smsq = 0
    smcs = 0
    for i in range(len(x)):
        smsq += x[i] ** 2
        smcs += np.cos(2 * np.pi * x[i])

    inx = 1 / len(x)
    return -20 * np.exp(-0.2 * np.sqrt(inx * smsq)) - np.exp(inx * smcs) + 20 + np.e


def rastrgin(x):
    sum = 0
    for i in range(len(x)):
        smsq = x[i] ** 2
        smcs = np.cos(2 * np.pi * x[i])
        sum = sum + smsq - 10 * smcs + 10
    return sum


def expanded_schaffers_f6(x):
    sum = 0
    for i in range(len(x)):
        if i != len(x) - 1:
            t = x[i] ** 2 + x[i + 1] ** 2
            sum += 0.5 + (((np.sin(np.sqrt(t))) ** 2 - 0.5) / ((1 + 0.001 * t) ** 2))
        elif i == len(x) - 1:
            t = x[i] ** 2 + x[0] ** 2
            sum += 0.5 + (((np.sin(np.sqrt(t))) ** 2 - 0.5) / ((1 + 0.001 * t) ** 2))
    return sum


def modified_schwefel(x):
    sm = 0
    for i in range(len(x)):
        z = x[i] + 420.9687462275036
        if np.abs(z) <= 500:
            sm += z * np.sin(np.sqrt(np.abs(z)))
        elif z > 500:
            sm += (500 - np.mod(z, 500)) * np.sin(np.sqrt(np.abs((500 - np.mod(z, 500))))) - ((z - 500) ** 2) / (
                        10000 * len(x))
        elif z < -500:
            sm += (np.mod(np.abs(z), 500) - 500) * np.sin(np.sqrt(np.abs((np.mod(np.abs(z), 500) - 500)))) - (
                    (z + 500) ** 2) / (10000 * len(x))
    return 418.9829 * len(x) - sm




def weierstrass(x):
    kcs = 0
    ksm = 0
    sm = 0
    for i in range(len(x)):
        for k in range(21):
            ak = 0.5 ** k
            bk = 2 * np.pi * (3 ** k)
            kcs += ak * np.cos((x[i] + 0.5) * bk)
            ksm += ak * np.cos((0.5 * bk))
        sm += kcs
    return sm - len(x) * ksm


def griewank(x):
    sm = 0
    cs = 1
    for i in range(len(x)):
        factor = 1 / 4000
        sm += factor * (x[i] ** 2)
        cs *= np.cos(x[i] / np.square(i + 1))
    return sm - cs + 1



def happy_cat(x):
    sm = 0
    smm = 0
    for i in range(len(x)):
        sm += (x[i]) ** 2
        smm += x[i]
    return (np.abs(sm - len(x))) ** 0.25 + (0.5 * sm + smm) / len(x) + 0.5



def func(n, X, op):
    if n == 1:
        if op == 'initial':
            return [-100, 100]
        else:
            return bent_cigar(X)

    elif n == 3:
        if op == 'initial':
            return [-100, 100]
        else:
            return zakharov(X)

    elif n == 4:
        if op == 'initial':
            return [-100, 100]
        else:
            return rosenbrock(X)

    elif n == 5:
        if op == 'initial':
            return [-100, 100]
        else:
            return rastrgin(X)

    elif n == 6:
        if op == 'initial':
            return [-100, 100]
        else:
            return expanded_schaffers_f6(X)

    elif n == 7:
        if op == 'initial':
            return [-100, 100]
        else:
            return 0

    elif n == 8:
        if op == 'initial':
            return [-100, 100]
        else:
            return 0

    elif n == 9:
        if op == 'initial':
            return [-100, 100]
        else:
            return levy(X)

    elif n == 10:
        if op == 'initial':
            return [-100, 100]
        else:
            return modified_schwefel(X)

    elif n == 11:
        if op == 'initial':
            return [-100, 100]
        else:
            return high_conditioned_elliptic(X)


    elif n == 12:
        if op == 'initial':
            return [-100, 100]
        else:
            return discus(X)

    elif n == 13:
        if op == 'initial':
            return [-100, 100]
        else:
            return ackley(X)


    elif n == 14:
        if op == 'initial':
            return [-100, 100]
        else:
            return weierstrass(X)


    elif n == 15:
        if op == 'initial':
            return [-100, 100]
        else:
            return griewank(X)

    elif n == 16:
        if op == 'initial':
            return [-100, 100]
        else:
            return katsuura(X)

    elif n == 17:
        if op == 'initial':
            return [-100, 100]
        else:
            return happy_cat(X)

    elif n == 18:
        if op == 'initial':
            return [-100, 100]
        else:
            return h_g_bat(X)

    elif n == 19:
        if op == 'initial':
            return [-100, 100]
        else:
            return 0

    elif n == 20:
        if op == 'initial':
            return [-100, 100]
        else:
            return schaffers_f7

    elif n == 1000:
        if op == 'initial':
            return [-100, 100]
        else:
            sum1 = 0
            for i in range(len(X)):
                sum1 += 0.5 * X[i]
            sum = 0
            for i in range(len(X)):
                sum += X[i] ** 2
            # f2
            summry = sum + sum1 ** 2 + sum1 ** 4

            sum2 = 0
            # f4
            for i in range(len(X)):
                sum2 += X[i] ** 2 - 10 * np.cos(2 * np.pi * X[i]) + 10

            sum3 = 0
            # f3
            for i in range(len(X) - 1):
                sum3 += 100 * (X[i + 1] - X[i] ** 2) ** 2 + (X[i] - 1) ** 2

            return 0.2 * summry + 0.4 * sum2 + 0.4 * sum3





    else:
        return