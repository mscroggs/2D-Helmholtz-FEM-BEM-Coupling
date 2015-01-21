from config import worldwide
from core import solve
from core import boundary
import cmath
from time import time
from math import log
from os.path import join
from numpy import array,linspace
from scipy import linalg
import matplotlib.pylab as plt

t00=time()

k=3
kc=k.conjugate()

dom=boundary.Domain([0,20])

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

xvgraph=[]
yvgraph=[]
yv2graph=[]

xvgraph=[a*10 for a in range(1,101)]

for igraph in xvgraph:
    err,timer=solve.helmholtz(int(igraph),f,dom,k,actual,img_show=False)

    yvgraph.append(err)
    yv2graph.append(timer)

t11=time()
print("Total time taken: %f seconds" %(t11-t00))

plt.clf()
plt.plot(xvgraph,yvgraph)
plt.title("Error")
plt.xlabel("Number of FEM nodes")
plt.ylabel("Max error")
plt.savefig(join(worldwide.PATH,'graph-error-inf.png'))
plt.show()

logxv=[log(a) for a in xvgraph]
logyv=[log(a) for a in yvgraph]

plt.clf()
plt.plot(logxv,logyv)
plt.title("Error")
plt.xlabel("log(Number of FEM nodes)")
plt.ylabel("log(Max error)")
plt.savefig(join(worldwide.PATH,'graph-error-inf-log.png'))
plt.show()

plt.clf()
plt.plot(xvgraph,yv2graph)
plt.title("Time Taken")
plt.xlabel("Number of FEM Nodes")
plt.ylabel("Time taken (seconds)")
plt.savefig(join(worldwide.PATH,'graph-time-taken-inf.png'))
plt.show()

