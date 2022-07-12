from tools.functions import Refractive_Index, SphericalMie, a_Lb_L, crosss
from tools.misc import Input

def Calculate_Individual():

    radius, PType, PFiles, MType, MFiles, eVmin, eVmax, wl_step = Input()
    n_particle, k_particle = Refractive_Index(PFiles, eVmin, eVmax, wl_step)
    n_matrix, k_matrix = Refractive_Index(MFiles, eVmin, eVmax, wl_step)
    a_L, b_L, L = a_Lb_L(n_particle, n_matrix, radius)
    hei = crosss(a_L, b_L, L, n_particle)
