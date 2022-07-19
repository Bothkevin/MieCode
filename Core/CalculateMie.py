from tools.functions import Refractive_Index, SphericalMie, a_Lb_L, crosss
from tools.misc import Input
from tools.functions_core_shell import a_Lb_L_core_shell

def Calculate_Individual():

    radius, PType, PFiles, MType, MFiles, eVmin, eVmax, wl_step = Input()
    n_particle, k_particle = Refractive_Index(PFiles, eVmin, eVmax, wl_step)
    n_matrix, k_matrix = Refractive_Index(MFiles, eVmin, eVmax, wl_step)
    n_particle[:,0,0] = n_particle[:,0,0]*1E-9
    k_particle[:,0,0] = k_particle[:,0,0]*1E-9
    n_matrix[:,0,0] = n_matrix[:,0,0]*1E-9
    k_matrix[:,0,0] = k_matrix[:,0,0]*1E-9
    a_L, b_L, L= a_Lb_L_core_shell(n_particle, k_particle, n_matrix, k_matrix, radius)
    #a_L, b_L, L= a_Lb_L(n_particle, k_particle, n_matrix,radius)
    hei = crosss(a_L, b_L, L, n_particle, n_matrix ,radius)
