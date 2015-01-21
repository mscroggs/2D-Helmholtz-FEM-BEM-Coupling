def helmholtz(phi_i,phi_j,dom,k):
    return -1j*k*(phi_i(dom[0])*phi_j(dom[0]).conjugate())-1j*k*(phi_i(dom[1])*phi_j(dom[1]).conjugate())
