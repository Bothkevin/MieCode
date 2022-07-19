import pandas as pd
import numpy as np
import scipy.constants as spc
from tools.misc import extrapolate, interpolate
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt

def Loadfile(file):

    # Returns the n and k of a csv file structured like the ones from
    # refractiveindex.info. wavelengths will be in nm

    data = pd.read_csv(file)

    k_start = data[data['n'] == 'k']

    if k_start.index > 1:
        ndata = pd.DataFrame(data[0:k_start.index[0]]).to_numpy(dtype = float)
        kdata = pd.DataFrame(data[k_start.index[0]+1:]).to_numpy(dtype = float)
    else:
        ndata = data.to_numpy(dtype = float)
        kdata = np.zeros(ndata.shape)
        kdata[:,0] = ndata[:,0]

    ndata[:,0] = ndata[:,0] * 1E3
    kdata[:,0] = kdata[:,0] * 1E3

    return ndata, kdata

def Refractive_Index(Files, eVmin, eVmax, wl_step):

    wlmin = int((spc.c * spc.h)/(eVmax * spc.e) * 1E9 - 1)
    wlmax = int((spc.c * spc.h)/(eVmin * spc.e) * 1E9 + 1)

    steps = int((wlmax - wlmin)/wl_step)

    ndata = np.zeros([steps,2,len(Files)])
    kdata = np.zeros([steps,2,len(Files)])
    filecounter = 0
    for file in Files:

        ndatatmp, kdatatmp = Loadfile(file)

        if ndatatmp[-1,0] < wlmax or ndatatmp[0,0] > wlmin:
            ntmp = extrapolate(ndatatmp, wlmin, wlmax)
            ktmp = extrapolate(kdatatmp, wlmin, wlmax)

            #This is for quality control at the moment TODO: remove
            plt.plot(ndatatmp[:,0], ndatatmp[:,1],color="r",linestyle='dashed',label = 'n input')
            plt.plot(ntmp[:,0], ntmp[:,1],color="g",linestyle='dashed',label = 'n extra')
            plt.legend()
            plt.savefig('Output\\n extrapolation of ' + file.split('\\')[-1].split('.')[0]+'.png', dpi = 300)
            plt.close()
            plt.plot(kdatatmp[:,0], kdatatmp[:,1],color="r",linestyle='dashed',label = 'k input')
            plt.plot(ktmp[:,0], ktmp[:,1],color="g",linestyle='dashed',label = 'k extra')
            plt.legend()
            plt.savefig('Output\\k extrapolation of ' + file.split('\\')[-1].split('.')[0]+'.png', dpi = 300)
            plt.close()

            ndata[:,:,filecounter] = interpolate(ntmp, wlmin, wlmax, wl_step)
            kdata[:,:,filecounter] = interpolate(ktmp, wlmin, wlmax, wl_step)

        else:
            ndata[:,:,filecounter] = interpolate(ndatatmp, wlmin, wlmax, wl_step)
            kdata[:,:,filecounter] = interpolate(kdatatmp, wlmin, wlmax, wl_step)

        #This is for quality control at the moment TODO: remove
        plt.plot(ndatatmp[:,0], ndatatmp[:,1],color="r",linestyle='dashed',label = 'n input')
        plt.plot(ndata[:,0,filecounter], ndata[:,1,filecounter],color="g",linestyle='dashed',label = 'n inter')
        plt.legend()
        plt.savefig('Output\\n interpolation of ' + file.split('\\')[-1].split('.')[0]+'.png', dpi = 300)
        plt.close()
        plt.plot(kdatatmp[:,0], kdatatmp[:,1],color="r",linestyle='dashed',label = 'k input')
        plt.plot(kdata[:,0,filecounter], kdata[:,1,filecounter],color="g",linestyle='dashed',label = 'k inter')
        plt.legend()
        plt.savefig('Output\\k interpolation of ' + file.split('\\')[-1].split('.')[0]+'.png', dpi = 300)
        plt.close()
    return ndata, kdata

def SphericalMie(z, L):
    spher_jn = spherical_jn(L,z,derivative=False)
    spherder_jn = spherical_jn(L,z,derivative=True)
    spher_yn = spherical_yn(L,z,derivative=False)
    spherder_yn = spherical_yn(L,z,derivative=True)

    riccarr_jn_data = z * spher_jn
    riccarrder_jn_data = spher_jn + z * spherder_jn
    riccarr_yn_data = z * spher_yn
    riccarrder_yn_data =  spher_yn + z * spherder_yn

    return riccarr_jn_data, riccarrder_jn_data, riccarr_yn_data, riccarrder_yn_data

def a_Lb_L(n, k, n_m, k_m,R):
    q = False
    filecounter = 0
    L = 1
    n_complex = n + 1j*k
    m =n_m[:, 1, filecounter]/ n_complex[:, 1,filecounter]

    x = np.abs(1/n[:,0,filecounter])*R
    a_Llist = []
    b_Llist = []
    while q == False:
        jn_Lmx = SphericalMie(m*x,L)[0]
        jnder_Lmx = SphericalMie(m*x,L )[1]
        jn_Lx = SphericalMie(x,L)[0]
        jnder_Lx = SphericalMie(x,L)[1]
        yn_Lx = SphericalMie(x,L)[2]
        ynder_Lx = SphericalMie(x,L)[3]

        nn_Lx = jn_Lx - 1j*yn_Lx
        nnder_Lx = jnder_Lx - 1j*ynder_Lx

        a_L = (m*jn_Lmx*jnder_Lx-jnder_Lmx*jn_Lx)/(m*jn_Lmx*nnder_Lx-jnder_Lmx*nn_Lx)
        b_L = (jn_Lmx*jnder_Lx-m*jnder_Lmx*jn_Lx)/(jn_Lmx*nnder_Lx-m*jnder_Lmx*nn_Lx)
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


    return a_L, b_L, L


def crosss(a_L, b_L, L, n, n_m, Radius):
    sigma_ext = 0
    sigma_sca = 0
    #print(Radius)
    #print(n[:,0, 0])
    for i in range(1,(L)):
        sigma_ext += (2*(i)+1)*(np.real(a_L[i-1]+b_L[i-1]))
        sigma_sca += (2*(i)+1)*(np.conjugate(a_L[i-1])*(a_L[i-1])+np.conjugate(b_L[i-1])*(b_L[i-1]))
    #print(f' ext {sigma_ext}')
    geo_cross = np.pi*Radius**2
    print(f'radius {Radius}')
    #print(k)
    kvec = np.abs((2*np.pi*n_m[:,1,0])/n_m[:,0,0])
    plt.figure(figsize = (15,5))
    plt.plot((spc.h*spc.c)/(n[:,0, 0]*spc.e),-2*np.pi/((kvec)**2*geo_cross)*sigma_ext,color="hotpink",linestyle='solid',label = 'extinction')
    plt.plot((spc.h*spc.c)/(n[:,0, 0]*spc.e),2*np.pi/((kvec)**2*geo_cross)*np.real(sigma_sca),color="plum",linestyle='dashed',label = 'scattering')
    plt.plot((spc.h*spc.c)/(n[:,0, 0]*spc.e),2*np.pi/((kvec)**2*geo_cross)*(sigma_ext - np.real(sigma_sca)) ,color="cornflowerblue",linestyle='dotted',label = 'absorption')
    plt.xlabel(f'$Energi [eV]$')
    plt.ylabel('$ \dfrac{ \sigma }{ \sigma_{geo} }$')
    plt.legend()
    plt.show()

    return print('hei')
