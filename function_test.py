import math
import numpy as np


def bent_cigar(X):
    # f1
    sum = X[0] * X[0]
    for i in range(1,len(X)):
        sum = sum + (X[i]**2)*10e6
    return sum

def zakharov(X):
    # f3
    sum1 = 0
    for i in range(len(X)):
        sum1 += 0.5 * X[i]
    sum = 0
    for i in range(len(X)):
        sum += X[i] ** 2
    return sum + sum1 ** 2 + sum1 ** 4

def rosenbrock(X):
    # f4
    sum = 0
    for i in range(len(X)-1):
        sum = sum + 100*(X[i+1] - X[i]**2)**2 + (X[i]-1)**2
    return sum

def rastrgin(x) :
    # f5
    sum = 0
    for i in range(len(x)):
        smsq = x[i] ** 2
        smcs = np.cos(2*np.pi * x[i])
        sum = sum + smsq - 10*smcs + 10
    return sum

def expanded_schaffers_f6(x) :
    # f6
    sum = 0
    for i in range(len(x)):
        if i != len(x)-1 :
            t = x[i] ** 2 + x[i+1] ** 2
            sum += 0.5 + (((np.sin(np.sqrt(t)))** 2-0.5)/((1 + 0.001*t)** 2))
        elif i == len(x)-1:
            t = x[i] ** 2 + x[0] ** 2
            sum += 0.5 + (((np.sin(np.sqrt(t)))** 2-0.5)/((1 + 0.001*t)** 2))
    return sum

def levy(x):
    # f9
    sum = 0
    for i in range(len(x)-1):
        w = 1 + (x[i]-1)/4
        sum +=((w-1)**2)*(1+10*((np.sin(np.pi*w+1))**2))
        # (np.sin(np.pi*(1 + (x[0]-1)/4)))**2 + sum + ((1 + (x[-1]/4) - 1)**2)*(1+np.sin(2*np.pi*(1 + (x[-1]-1)/4))**2)
    return (np.sin(np.pi*(1 + (x[0]-1)/4)))**2 + sum + ((1 + (x[i]/4) - 1)**2)*(1+np.sin(2*np.pi*(1 + (x[i]-1)/4))**2)

def modified_schwefel(x):
    # f10
    sm = 0
    for i in range(len(x)):
        z = x[i]+420.9687462275036
        if np.abs(z) <= 500:
            sm += z * np.sin(np.sqrt(np.abs(z)))
        elif z > 500:
            sm += (500 - np.mod(z, 500))*np.sin(np.sqrt(np.abs((500 - np.mod(z, 500)))))-((z-500)**2)/(10000*len(x))
        elif z < -500:
            sm += (np.mod(np.abs(z), 500) - 500) * np.sin(np.sqrt(np.abs(( np.mod(np.abs(z), 500)- 500)))) - (
                        (z + 500) ** 2) / (10000 * len(x))
    return 418.9829*len(x)-sm

def high_conditioned_elliptic(x):
    # f11
    sm = 0
    for i in range(len(x)):
          sm += (x[i]**2)*1000000**(i/(len(x)-1))
    return sm

def discus(x):
    # f12
    sm = 0
    for i in range(1,len(x)):
          sm += x[i]**2
    return sm+(x[0]**2)*1000000

def ackley(x):
    # f13
    smsq = 0
    smcs = 0
    for i in range(len(x)):
        smsq += x[i] ** 2
        smcs += np.cos(2*np.pi * x[i])
    
    inx = 1/len(x)
    return -20*np.exp(-0.2*np.sqrt(inx*smsq)) - np.exp(inx*smcs) + 20 + np.e


def weierstrass1(x):
    kcs=0
    ksm=0
    sm=0
    for i in range(len(x)):

        for k in range(21):
            ak = 0.5 ** k
            bk = 2*np.pi * (3 ** k)
            kcs += ak * np.cos((x[i] + 0.5) * bk)
            ksm += ak * np.cos((0.5 * bk))
        sm +=kcs
    return sm - len(x)*ksm

def weierstrass(x):
    # f14
    sum1 = np.zeros_like(x)
    sum2 = 0
    for k in range(20):
        sum1 = sum1 + (0.5**(k+1))*np.cos(2*np.pi*(x+0.5)*(3**(k+1)))
        sum2 = sum2 + (0.5**(k+1))*np.cos(2*np.pi*(3**(k+1))*0.5)
    res = np.sum(sum1)-len(x)*sum2
    return res


def griewank(x):
    # f15
    sm = 0
    cs = 1
    for i in range(len(x)):
        factor = 1 / 4000
        sm += factor*(x[i]**2)
        cs *= np.cos(x[i]/np.square(i+1))
    return sm-cs+1

def katsuura(x):
    # f16
    cs = 1
    D = len(x)
    d = 10/(D**2)
    sum1 = np.zeros_like(x)
    for j in range(1,33):
        sum1 += (np.abs((2**j)*x-np.round((2**j)*x)) / (2**j))
    li = list(range(1,D+1))
    sum1 = (sum1 * li + 1)**(10/(D**1.2))
    cs = np.cumprod(sum1)[-1]
    return (d*cs)- d

def happy_cat(x):
    # f17
    sm = 0
    smm = 0
    for i in range(len(x)):
        sm += (x[i])**2
        smm += x[i]
    return (np.abs(sm-len(x)))**0.25 + (0.5 * sm + smm)/len(x) + 0.5

def h_g_bat(x):
    # f18
    sm = 0
    smm = 0
    for i in range(len(x)):
        sm += (x[i]) ** 2
        smm += x[i]
    return (np.abs(sm**2-sm))**0.5 + (0.5 * sm + smm)/len(x) + 0.5

def schaffers_f7(x):
    # f20
    sm = 0
    for i in range(len(x)-1):
        s = np.square((x[i])**2+(x[i+1])**2)
        sm += np.square(s)*(np.sin(50*(s)**0.2))+1
    return (sm/(len(x)-1))**2

def Hybrid_Function_1(x):
    # f21
    sum1 = zakharov(x)
    sum2 = rosenbrock(x)
    sum3 = rastrgin(x)
    return 0.2*sum1 + 0.4*sum2 +0.4*sum3

def Hybrid_Function_2(x):
    # f22
    sum1 = high_conditioned_elliptic(x)
    sum2 = modified_schwefel(x)
    sum3 = bent_cigar(x)
    return 0.3*sum1 + 0.3*sum2 +0.4*sum3

def Hybrid_Function_3(x):#缺f7
    # f23
    # sum1 = high_conditioned_elliptic(x)
    # sum2 = discus(x)
    # sum3 = bent_cigar(x)
    return 0

def Hybrid_Function_4(x):
    # f24
    sum1 = high_conditioned_elliptic(x)
    sum2 = ackley(x)
    sum3 = schaffers_f7(x)
    sum4 = rastrgin(x)
    return 0.2*sum1 + 0.2*sum2 +0.2*sum3+0.4*sum4

def Hybrid_Function_5(x):
    # f25
    sum1 = bent_cigar(x)
    sum2 = h_g_bat(x)
    sum3 = rastrgin(x)
    sum4 = rosenbrock(x)
    return 0.2*sum1 + 0.2*sum2 +0.3*sum3+0.3*sum4

def Hybrid_Function_6(x):
    # f26
    sum1 = expanded_schaffers_f6(x)
    sum2 = h_g_bat(x)
    sum3 = rosenbrock(x)
    sum4 = modified_schwefel(x)
    return 0.2*sum1 + 0.2*sum2 +0.3*sum3+0.3*sum4

def Hybrid_Function_7(x):#缺f19
    # f27
    sum1 = katsuura(x)
    sum2 = ackley(x)
    sum3 = expanded_schaffers_f6(x)
    sum4 = modified_schwefel(x)
    sum5 = rastrgin(x)
    return 0.1*sum1 + 0.2*sum2 +0.2*sum3+0.2*sum4+0.3*sum5

def Hybrid_Function_8(x):
    # f28
    sum1 = high_conditioned_elliptic(x)
    sum2 = ackley(x)
    sum3 = rastrgin(x)
    sum4 = h_g_bat(x)
    sum5 = discus(x)
    return 0.2*sum1 + 0.2*sum2 +0.2*sum3+0.2*sum4+0.2*sum5

def Hybrid_Function_9(x):#缺f19
    # f29
    sum1 = bent_cigar(x)
    sum2 = rastrgin(x)
    sum3 = expanded_schaffers_f6(x)
    sum4 = weierstrass(x)
    sum5 = expanded_schaffers_f6(x)
    return 0.2*sum1 + 0.2*sum2 +0.2*sum3+0.2*sum4+0.2*sum5

def Hybrid_Function_10(x):
    # f30
    sum1 = happy_cat(x)
    sum2 = katsuura(x)
    sum3 = ackley(x)
    sum4 = rastrgin(x)
    sum5 = modified_schwefel(x)
    sum6 = schaffers_f7(x)
    return 0.1*sum1 + 0.1*sum2 +0.2*sum3+0.2*sum4+0.2*sum5+0.2*sum6


def func(n,X,op):
    if n==1:
        if op == 'initial':
            return [-100,100]
        else:
            return bent_cigar(X)

    elif n==3:
        if op == 'initial':
            return [-100,100]
        else:
            return zakharov(X)

    elif n==4:
        if op == 'initial':
            return [-100,100]
        else:
            return rosenbrock(X)

    elif n==5:
        if op == 'initial':
            return [-5,5]
        else:
            return rastrgin(X)

    elif n==6:
        if op == 'initial':
            return [-100,100]
        else:
            return expanded_schaffers_f6(X)
    
    elif n==7:
        if op == 'initial':
            return [-100,100]
        else:
            return 0
    
    elif n==8:
        if op == 'initial':
            return [-100,100]
        else:
            return 0

    elif n==9:
        if op == 'initial':
            return [-100,100]
        else:
            return levy(X)

    elif n==10:
        if op == 'initial':
            return [-100,100]
        else:
            return modified_schwefel(X)  

    elif n==11:
        if op == 'initial':
            return [-100,100]
        else:
            return high_conditioned_elliptic(X)   

    elif n==12:
        if op == 'initial':
            return [-100,100]
        else:
            return discus(X)

    elif n==13:
        if op == 'initial':
            return [-32,32]
        else:
            return ackley(X)

    elif n==14:
        if op == 'initial':
            return [-100,100]
        else:
            return weierstrass(X)

    elif n==15:
        if op == 'initial':
            return [-600,600]
        else:
            return griewank(X)

    elif n==16:
        if op == 'initial':
            return [-100,100]
        else:
            return katsuura(X)

    elif n==17:
        if op == 'initial':
            return [-100,100]
        else:
            return happy_cat(X)

    elif n==18:
        if op == 'initial':
            return [-100,100]
        else:
            return h_g_bat(X)

    elif n==19:
        if op == 'initial':
            return [-100,100]
        else:
            return 0

    elif n==20:
        if op == 'initial':
            return [-100,100]
        else:
            return schaffers_f7(X)

    elif n==21:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_1(X)

    elif n==22:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_2(X)

    elif n==23:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_3(X)
    
    elif n==24:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_4(X)

    elif n==25:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_5(X)

    elif n==26:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_6(X)

    elif n==27:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_7(X)

    elif n==28:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_8(X)

    elif n==29:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_9(X)

    elif n==30:
        if op == 'initial':
            return [-100,100]
        else:
            return Hybrid_Function_10(X)



    elif n==1000:
        if op == 'initial':
            return [-100,100]
        else:
            sum1 = 0
            for i in range(len(X)):
                sum1 += 0.5 * X[i]
            sum = 0
            for i in range(len(X)):
                sum += X[i] ** 2
            #f2
            summry = sum + sum1 ** 2 + sum1 ** 4

            sum2 = 0
            #f4
            for i in range(len(X)):
                sum2 += X[i] ** 2 - 10 * np.cos(2 * np.pi * X[i]) + 10
            
            sum3 = 0
            #f3
            for i in range(len(X) - 1):
                sum3 += 100 * (X[i+1] - X[i] ** 2)**2 + (X[i] - 1)**2
            
            return 0.2 * summry + 0.4 * sum2 + 0.4 * sum3
    




    else:
        return