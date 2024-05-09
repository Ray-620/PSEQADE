import math
import numpy as np
import cec2017.functions as functions
import cec2017.basic as basic


def func(n,x,op):
    X = [x]

    if n==1:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f1
            return f(X)[0]
    elif n==2:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f2
            return f(X)[0]

    elif n==3:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f3
            return f(X)[0]

    elif n==4:
        if op == 'initial':
            return [-10,10]
        else:
            f = functions.f4
            return f(X)[0]
    
    elif n==5:
        if op == 'initial':
            return [-10,10]
        else:
            f = functions.f5
            return f(X)[0]
    
    elif n==6:
        if op == 'initial':
            return [-20,20]
        else:
            f = functions.f6
            return f(X)[0]
    
    elif n==7:
        if op == 'initial':
            return [-50,50]
        else:
            f = functions.f7
            return f(X)[0]

    elif n==8:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f8
            return f(X)[0]
        
    elif n==9:
        if op == 'initial':
            return [-10,10]
        else:
            f = functions.f9
            return f(X)[0]


    elif n==10:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f10
            return f(X)[0]

    elif n==11:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f11
            return f(X)[0]

    elif n==12:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f12
            return f(X)[0]

    elif n==13:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f13
            return f(X)[0]

    elif n==14:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f14
            return f(X)[0]

    elif n==15:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f15
            return f(X)[0]

    elif n==16:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f16
            return f(X)[0]

    elif n==17:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f17
            return f(X)[0]

    elif n==18:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f18
            return f(X)[0]

    elif n==19:
        if op == 'initial':
            return [-50,50]
        else:
            f = functions.f19
            return f(X)[0]

    elif n==20:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f20
            return f(X)[0]

    elif n==21:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f21
            return f(X)[0]

    elif n==22:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f22
            return f(X)[0]

    elif n==23:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f23
            return f(X)[0]
    
    elif n==24:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f24
            return f(X)[0]

    elif n==25:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f25
            return f(X)[0]

    elif n==26:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f26
            return f(X)[0]

    elif n==27:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f27
            return f(X)[0]
    
    elif n==28:
        if op == 'initial':
            return [-50,50]
        else:
            f = functions.f28
            return f(X)[0]
        
    elif n==29:
        if op == 'initial':
            return [-100,100]
        else:
            f = functions.f29
            return f(X)[0]

    elif n==30:
        if op == 'initial':
            return [-100,100]
        else:
            f = basic.zakharov
            return f(X)[0]
    

    else:
        return