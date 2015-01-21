from config import worldwide
from core import boundary
from core import solve
import cmath
from numpy import array
from scipy import linalg
#import matplotlib.pylab as plt

k=3
kc=k.conjugate()
dom=boundary.Domain([0,20])

#find the actual solution
A=array([[0]*4]*4,dtype=complex)
b=array([0]*(4),dtype=complex)
#value at dom[0]
A[0,0]=cmath.exp(1j*kc*dom[0])
A[0,2]=-cmath.exp(1j*kc*dom[0])
A[0,3]=-cmath.exp(-1j*kc*dom[0])
b[0]=float(dom[0])/kc**2
#derivative at dom[0]
A[1,0]=1j*kc*cmath.exp(1j*kc*dom[0])
A[1,2]=-1j*kc*cmath.exp(1j*kc*dom[0])
A[1,3]=1j*kc*cmath.exp(-1j*kc*dom[0])
b[1]=1.0/kc**2
#value at dom[1]
A[2,1]=cmath.exp(-1j*kc*dom[1])
A[2,2]=-cmath.exp(1j*kc*dom[1])
A[2,3]=-cmath.exp(-1j*kc*dom[1])
b[2]=float(dom[1])/kc**2
#derivative at dom[1]
A[3,1]=-1j*kc*cmath.exp(-1j*kc*dom[1])
A[3,2]=-1j*kc*cmath.exp(1j*kc*dom[1])
A[3,3]=1j*kc*cmath.exp(-1j*kc*dom[1])
b[3]=1.0/kc**2

coeffs=linalg.solve(A,b)

def actual(x):
    if x<dom[0]:
        return coeffs[0]*cmath.exp(1j*kc*x)
    elif x<dom[1]:
        return coeffs[2]*cmath.exp(1j*kc*x)+coeffs[3]*cmath.exp(-1j*kc*x)+x/kc**2
    else:
        return coeffs[1]*cmath.exp(-1j*kc*x)

def f(x):
    return -x

#inputs here are (number of FEM nodes,f,domain,k,actual solution)
solve.helmholtz(8,f,dom,k,actual)
