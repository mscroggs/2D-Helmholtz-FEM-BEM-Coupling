from core import boundary
from core import solve
import cmath
from numpy import array
from scipy import linalg
#import matplotlib.pylab as plt

k=11

def f(x):
    return -1

dom=boundary.Domain([0,20])

solve.helmholtz(10,f,dom,k,fix="AUTO")

