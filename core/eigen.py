from config import worldwide
from core import element
from core import boundary
import cmath
from core import DtN
from scipy import integrate
from scipy.misc import derivative
from scipy.linalg import eig
from numpy import array
from os.path import join

with open(join(worldwide.PATH,'eigen'),'w+') as f:
    pass

def fwrite(txt):
    with open(join(worldwide.PATH,'eigen'),'a+') as f:
        f.write(txt+"\n")

def helmholtz(Ninput,dom):
    def DtNmap(phi_i,phi_j):
        return DtN.helmholtz(phi_i,phi_j,dom,1)
    solve(Ninput,dom,DtNmap)

def solve(Ninput,dom,DtNmap):
    #prepare domain
    dom.extend(Ninput)

    #create blank matrices
    #A=[[K,-D],[O,I]]
    #B=[[O,M],[I,O]]
    A=array([[0]*(2*Ninput+4)]*(2*Ninput+4),dtype=complex)
    B=array([[0]*(2*Ninput+4)]*(2*Ninput+4),dtype=complex)

    #I
    for i in range(0,Ninput+2):
        A[Ninput+2+i,Ninput+2+i]=1
        B[Ninput+2+i,i]=1
    #FEM
    for i in range(0,Ninput+2):
        phi_i,dphi_i=element.new_linear(dom(i-1),dom(i),dom(i+1))
        int_over_me=[dom(i-1),dom(i+1)]
        xv=[]
        yv=[]
        yv2=[]

        def integrand(x): return phi_i(x)*f(x)

        for j in range(0,Ninput+2):
          if j>=i-1 and j<=i+1:
            Aentry=0
            phi_j,dphi_j=element.new_linear(dom(j-1),dom(j),dom(j+1))
            def integrand1(x): return dphi_i(x)*dphi_j(x)
            def integrand2(x): return phi_i(x)*phi_j(x)
            i1,error=integrate.quad(integrand1,int_over_me[0],int_over_me[1])
            i2,error=integrate.quad(integrand2,int_over_me[0],int_over_me[1])
            K=i1
            M=-i2
            D=DtNmap(phi_i,phi_j)

            #Amat=K-k**2*M-k*D

            #K
            A[i,j]=K
            #M
            B[i,Ninput+2+j]=M
            #D
            A[i,Ninput+2+j]=-D
    eigen=eig(A,B)
    for ei in eigen[0]:
        fwrite(str(ei))
