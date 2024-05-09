import math
import numpy as np
from cec2013.functions import CEC_functions
from cec2013lsgo.cec2013 import Benchmark


def func(n,X,op):
    bench = Benchmark()
    info = bench.get_info(n)

    if op == 'initial':
        return [info['lower'],info['upper']]
    else:
        fn = bench.get_function(n)
        return fn(X)


