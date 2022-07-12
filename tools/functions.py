import pandas as pd
import numpy as np
import scipy.constants as spc
from tools.misc import extrapolate, interpolate

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

def SphericalMie():
