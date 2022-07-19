import pandas as pd
import numpy as np
import scipy.constants as spc
from tools.misc import extrapolate, interpolate
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt
from tools.functions import SphericalMie

def a_Lb_L_core_shell(n, k, n_m, k_m, R):
    q = False
    filecounter = 0
    L = 1
    n_complex = n - 1j*k
    n_m_complex = n_m - 1j*k_m
    #print(f'n vanlig{n} n complex{n_complex}, k {k}')

    m = n_m_complex[:, 1, filecounter]/n_complex[:, 1,filecounter]
    kvec = (2*np.pi*n_complex[:,1,filecounter])/n[:,0,filecounter]

    #print(f'k-vector absoluttverdi kvadrat {np.sqrt(np.real(kvec*np.conjugate(kvec)))} ')
    #x = np.sqrt(np.real(kvec*np.conjugate(kvec)))*R
    x = kvec*R
    a_Llist = []
    b_Llist = []

    while q == False:
        jn_Lmx, jnder_Lmx, yn_Lmx, ynder_Lmx  = SphericalMie(m*x,L)
        jn_Lx, jnder_Lx, yn_Lx, ynder_Lx  = SphericalMie(x,L)

        xin_Lx = jn_Lx - 1j*yn_Lx
        xinder_Lx = jnder_Lx -  1j*ynder_Lx
        xin_Lmx = jn_Lmx - 1j*yn_Lmx
        xinder_Lmx = jnder_Lmx - 1j*ynder_Lmx

        etan_Lx = jn_Lx + 1j*yn_Lx
        etander_Lx = jnder_Lx + 1j*ynder_Lx
        etan_Lmx = jn_Lmx + 1j*yn_Lmx
        etander_Lmx = jnder_Lmx + 1j*ynder_Lmx

        T_L0 = 0
        S_L0 = 0

        #T_L1 = -(m*jn_Lmx*(jnder_Lx+T_L0*ynder_Lx)-jnder_Lmx*(jn_Lx+T_L0*yn_Lx))/(m*yn_Lmx*(jnder_Lx+T_L0*ynder_Lx)-ynder_Lmx*(jn_Lx+T_L0*yn_Lx))
        #S_L1 = -(jn_Lmx*(jnder_Lx+S_L0*ynder_Lx)-m*jnder_Lmx*(jn_Lx+S_L0*yn_Lx))/(yn_Lmx*(jnder_Lx+S_L0*ynder_Lx)-m*ynder_Lmx*(jn_Lx+S_L0*yn_Lx))
        #print(f'first T{np.shape(T_L1)}')

        h = 1
        T_L = np.zeros((h+1,len(x)),dtype=np.complex_)
        S_L = np.zeros((h+1,len(x)),dtype=np.complex_)
        #print(T_L)
        T_L[0,:] = T_L0
        S_L[0,:] = S_L0
        #print(T_L)
        for i in range(h):
            #print('hola hoho')
            #mer riktig
            T_L[i+1,:] = -(m*jn_Lmx*(jnder_Lx+T_L[i,:]*ynder_Lx)-jnder_Lmx*(jn_Lx+T_L[i,:]*yn_Lx))/(m*yn_Lmx*(jnder_Lx+T_L[i,:]*ynder_Lx)-ynder_Lmx*(jn_Lx+T_L[i,:]*yn_Lx))
            S_L[i+1,:] = -(jn_Lmx*(jnder_Lx+S_L[i,:]*ynder_Lx)-m*jnder_Lmx*(jn_Lx+S_L[i,:]*yn_Lx))/(yn_Lmx*(jnder_Lx+S_L[i,:]*ynder_Lx)-m*ynder_Lmx*(jn_Lx+S_L[i,:]*yn_Lx))


        #For for-loop
        #hvilket skall vi skal se på lol
        s = h
        a_L = -(m*jn_Lmx*(jnder_Lx+T_L[s]*ynder_Lx)-jnder_Lmx*(jn_Lx+T_L[s]*yn_Lx))/(m*xin_Lmx*(jnder_Lx+T_L[s]*ynder_Lx)-xinder_Lmx*(jn_Lx+T_L[s]*yn_Lx))
        b_L = -(jn_Lmx*(jnder_Lx+S_L[s]*ynder_Lx)-m*jnder_Lmx*(jn_Lx+S_L[s]*yn_Lx))/(xin_Lmx*(jnder_Lx+S_L[s]*ynder_Lx)-m*xinder_Lmx*(jn_Lx+S_L[s]*yn_Lx))



        if any(np.isnan(a_L)) == True:
            q = True
        elif any(np.isnan(b_L)) == True:
            q = True
        else:
            L += 1
            a_Llist.append(a_L)
            b_Llist.append(b_L)
        a_L = np.asarray(a_Llist)
        b_L = np.asarray(b_Llist)
        #print(f'T{h} hmm{T_L[h,:]}')
        #print(f'S{h} hmm{S_L[:,h]}')
    return a_L, b_L, L
