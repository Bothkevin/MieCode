import os
import glob
import numpy as np
from scipy.interpolate import interp1d

def inputfiles():

    cwd = os.getcwd()
    inputdir = cwd + '\\InputFiles\\'
    files = glob.glob(inputdir + '*.csv')

    return files


def files_manager(Particles, Matrices):

    ParticleFiles = []
    MatrixFiles = []

    if not Particles:
        print(r'No particle selected!')

    elif not Matrices:
        print(r'No matrix selected!')

    else:
        input_files = inputfiles()
        structure = []

        for file in input_files:
            name = file.split('\\')[-1]
            structure.append(name.split('_')[0])

        for Particle in Particles:
            index = structure.index(Particle)
            ParticleFiles.append(input_files[index])

        for Matrix in Matrices:
            index = structure.index(Matrix)
            MatrixFiles.append(input_files[index])

        if len(ParticleFiles) != len(Particles):
            print('File - Particles mismatch!')

        if len(MatrixFiles) != len(Matrices):
            print('File - Matrices mismatch!')

    return ParticleFiles, MatrixFiles

def extrapolate(data, L_boundary, U_boundary):

    step = 1
    points = 4

    if data[0,0] > L_boundary:
        minrange = np.zeros([points,2])
        for i in range(0, points, 1):
            minrange[i,:] = data[0,:]
            data = np.delete(data,0,0)
        min_val = np.arange(L_boundary, minrange[points-1,0], 1)

        data_min_f = np.polyfit(minrange[:,0], minrange[:,1],2)
        data_min_p = np.poly1d(data_min_f)
        data_min = data_min_p(min_val)

        min_arr = np.transpose(np.stack((min_val,data_min),axis=0))
        extrapolate_temp = np.concatenate((min_arr, data))
    else:
        extrapolate_temp = data

    if extrapolate_temp[-1,0] < U_boundary:
        maxrange = np.zeros([points,2])
        for j in range(0, points, 1):
            maxrange[j-1,:] = data[-1,:]
            data = np.delete(data,-1,0)
        max_val = np.arange(U_boundary, maxrange[0,0],1)

        data_max_f = np.polyfit(maxrange[:,0], maxrange[:,1],2)
        data_max_p = np.poly1d(data_max_f)
        data_max = data_max_p(max_val)

        max_arr = np.transpose(np.stack((max_val, data_max), axis = 0))
        extrapolate_data = np.concatenate((extrapolate_temp,max_arr))
    else:
        extrapolate_data = extrapolate_temp

    return extrapolate_data
    

def interpolate(data, L_boundary, U_boundary, step):

    xval = np.arange(L_boundary, U_boundary, step, dtype = int)
    Interpolation_F = interp1d(data[:, 0], data[:, 1], kind = 'cubic')
    yval_interpolate = Interpolation_F(xval)

    interpolated_data = np.transpose(np.stack((xval, yval_interpolate), axis = 0))

    return interpolated_data


def Input():

    eVmin = 1
    eVmax = 3
    wl_step = 1

    R = 20*1E-9
    Type = 'Sphere'
    Particles = ['Au']

    MatrixType = 'Uniform'
    Matrices = ['STO']

    Particlefiles, Matrixfiles = files_manager(Particles,Matrices)

    return R, Type, Particlefiles, MatrixType, Matrixfiles, eVmin, eVmax, wl_step
