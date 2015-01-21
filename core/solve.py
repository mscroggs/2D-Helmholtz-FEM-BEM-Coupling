from config import worldwide
from core import boundary
from core import element
from core import DtN
import cmath
from time import time
from scipy import integrate
from scipy.misc import derivative
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.optimize import minimize_scalar as minimise_scalar
from os.path import join
from numpy import array,real,linspace
import matplotlib.pylab as plt

#def poisson(Ninput,f,dom,k,actual=None,img_show=True):
#    def y(x):
#        return -abs(x)/2
#
#    def DtNmap(phi_i,phi_j):
#        return DtN.poisson(phi_i,phi_j,dom,k)
#    name="Poisson"
#    return solve(Ninput,f,dom,k,y,DtNmap,name,actual,img_show)

def helmholtz(Ninput,f,dom,k,actual=None,img_show=True,fix=None):
    def y(x):
        return cmath.exp(1j*k*x)

    def DtNmap(phi_i,phi_j):
        return DtN.helmholtz(phi_i,phi_j,dom,k)
    name="Helmholtz"
    return solve(Ninput,f,dom,k,y,DtNmap,name,actual,img_show,fix)

def jump(f,x,dx):
    return f(x-dx)-f(x+dx)

def solve(Ninput,f,dom,k,y,DtNmap,name,actual=None,img_show=True,fix=None):
    t0=time()

    # Erik's 1D fixing
    if fix=="AUTO":
        fixt=k*dom.width
        fix=(6*cmath.cos(fixt)-6+fixt**2*cmath.cos(fixt)+2*fixt**2)/(12*(1-cmath.cos(fixt))**2)

    # set up mesh
    dom.extend(Ninput)

    # set up solution vectors
    u=[None]*(Ninput+2)

    # create matrix
    A=lil_matrix((Ninput+2,Ninput+2),dtype=complex)
    b=array([0]*(Ninput+2),dtype=complex)

    #FEM
    for i in range(0,Ninput+2):
        phi_i,dphi_i=element.new_linear(dom(i-1),dom(i),dom(i+1))
        int_over_me=[dom(i-1),dom(i+1)]
        xv=[]
        yv=[]
        yv2=[]

        def integrand(x): return phi_i(x)*f(x)
        bnew,error=integrate.quad(integrand,int_over_me[0],int_over_me[1])

        for j in range(0,Ninput+2):
          if j>=i-2 and j<=i+2:
            Aentry=0
            phi_j,dphi_j=element.new_linear(dom(j-1),dom(j),dom(j+1))
            def integrand1(x): return dphi_i(x)*dphi_j(x)
            def integrand2(x): return phi_i(x)*phi_j(x)
            i1,error=integrate.quad(integrand1,int_over_me[0],int_over_me[1])
            i2,error=integrate.quad(integrand2,int_over_me[0],int_over_me[1])
            Aentry+=i1
            Aentry-=k**2*i2
            Aentry+=DtNmap(phi_i,phi_j)
            if fix!=None:
                for ki in range(max(1,i-1),min(Ninput,i+2)):
                    Aentry+=fix*jump(dphi_i,dom(ki),dom.width/2)*jump(dphi_j,dom(ki),dom.width/2)*dom.width
                Aentry+=fix*dom.width*(dphi_i(dom[1])-1j*k*phi_i(dom[1]))*(dphi_j(dom[1]).conjugate()+1j*k*phi_j(dom[1]).conjugate())
                Aentry-=fix*dom.width*(dphi_i(dom[0])-1j*k*phi_i(dom[0]))*(dphi_j(dom[0]).conjugate()+1j*k*phi_j(dom[0]).conjugate())
                #FIX THESE TWO LINES
            A[i,j]=Aentry
        b[i]=bnew

    #solve
    A=A.tocsr()
    b=array(b)
    soln=spsolve(A,b)
    def solution(x):
        if x<dom[0]:
            return soln[0]*y(dom[0]-x)
        for i in range(1,len(u)):
            if x<dom(i):
                return (soln[i]-soln[i-1])*(x-dom(i-1))/(dom(i)-dom(i-1))+soln[i-1]
        return soln[-1]*y(x-dom[1])
    t1=time()

    #draw some pretty pictures
    if img_show:
        uMesh=linspace(dom[0]-dom.FEM_width/2,dom[1]+dom.FEM_width/2,3000)
        uPlot=[solution(p) for p in uMesh]
        if actual!=None:
            yv=[actual(p) for p in uMesh]
            yer=[abs(real(a-b)) for a,b in zip(uPlot,yv)]
        plt.clf()
        plt.plot(uMesh,uPlot)
        labels=['FEM BEM Coupling']
        if actual!=None:
            plt.plot(uMesh,yv)
            labels.append('Actual Solution')
        plt.legend(labels, loc='upper left',prop={'size':8})
        plt.grid(True)
        plt.title(name+" FEM in ("+str(dom[0])+","+str(dom[1])+") with "+str(Ninput)+" mesh nodes")
        plt.savefig(join(worldwide.PATH,'graph.png'))
        plt.show()
        plt.clf()
        maxyer=0
        if actual!=None:
            plt.title("Error for "+name+" FEM in ("+str(dom[0])+","+str(dom[1])+") with "+str(Ninput)+" mesh nodes")
            plt.plot(uMesh,yer)
            plt.legend(['Error'],loc='upper left',prop={'size':8})
            plt.grid(True)
            plt.savefig(join(worldwide.PATH,'graph-error.png'))
            plt.show()
            maxyer=max(yer)
    elif actual!=None:
        def yer(x):
            return -abs(real(solution(x)-actual(x)))
        maxyer=-minimise_scalar(yer).fun
        t1=time()

    print("Time taken for "+str(int(Ninput))+" nodes: %f seconds" %(t1-t0))
    if actual!=None:
        return maxyer,t1-t0
    else:
        return 0,t1-t0
