"""
Created on Thu Feb  7 10:00:12 2019
Some mathematical functions for paleochrono.
@author: parrenif
"""
import numpy as np
import math as m
from numba import jit, njit, prange, guvectorize, float64

def interp_lin_aver(x_out, x_in, y_in):
    """Return a linear interpolation of a (x_in,y_in) series at x_out abscissas with averaging."""
    y_out = np.nan*np.zeros(np.size(x_out)-1)
    if x_out[0] < min(x_in):
        x_mod = np.concatenate((np.array([x_out[0]]), x_in))
        y_mod = np.concatenate((np.array([y_in[0]]), y_in))
    else:
        x_mod = x_in+0
        y_mod = y_in+0
    if x_out[-1] > max(x_in):
        x_mod = np.concatenate((x_mod, np.array([x_out[-1]])))
        y_mod = np.concatenate((y_mod, np.array([y_in[-1]])))
    for i in range(np.size(x_out)-1):
        x_loc = x_mod[np.where(np.logical_and(x_mod > x_out[i], x_mod < x_out[i+1]))]
        x_loc = np.concatenate((np.array([x_out[i]]), x_loc, np.array([x_out[i+1]])))
        y_loc = np.interp(x_loc, x_mod, y_mod)
        y_out[i] = np.sum((y_loc[1:]+y_loc[:-1])/2*(x_loc[1:]-x_loc[:-1]))/(x_out[i+1]-x_out[i])
    return y_out

def interp_stair_aver(x_out, x_in, y_in):
    """Return a staircase interpolation of a (x_in,y_in) series at x_out abscissas with averaging.
    """
    x_mod = x_in+0
    y_mod = y_in+0
    if x_out[0] < x_in[0]:
        x_mod = np.concatenate((np.array([x_out[0]]), x_mod))
        y_mod = np.concatenate((np.array([y_in[0]]), y_mod))
    if x_out[-1] > x_in[-1]:
        x_mod = np.concatenate((x_mod, np.array([x_out[-1]])))
        y_mod = np.concatenate((y_mod, np.array([y_in[-1]])))
    y_int = np.cumsum(np.concatenate((np.array([0]), y_mod[:-1]*(x_mod[1:]-x_mod[:-1]))))
#Maybe this is suboptimal since we compute twice g(xp[i]):
    y_out = (np.interp(x_out[1:], x_mod, y_int)-np.interp(x_out[:-1], x_mod, y_int))/\
            (x_out[1:]-x_out[:-1])
    return y_out


def gaussian(x_in):
    """Return the value of the gaussian function (no multiplicative constant)
    at a given x_in abscissa."""
    return np.exp(-x_in**2/2)

def grid(para):
    start = para['start']
    end = para['end']
    try:
        nb_steps = para['nb_steps']
    except KeyError:
        resolution = para['resolution']
        nb_steps = m.floor((end-start)/resolution)
        end = start + resolution * nb_steps
    if para['type'] == 'regular':
        eps = (end-start)/nb_steps/2
        grid = np.arange(start, end+eps, (end-start)/nb_steps)
    elif para['type'] == 'linear':
        ratio = para['ratio']
        if ratio == None:
            ratio = 2/(nb_steps+1)
        eps = (1.-ratio)/nb_steps
        grid = np.arange(ratio, 2.-ratio+eps, (2.-2*ratio)/(nb_steps-1))
        grid = grid * (end-start)/nb_steps
        grid = np.cumsum(np.concatenate((np.array([start]), grid)))
    else:
        print('Type of grid not recognized.')
    try:
        inverted = para['inverted']
    except KeyError:
        inverted = False
    if inverted:
        grid = grid[::-1]
        grid = grid[:-1]-grid[1:]
        grid = np.cumsum(np.concatenate((np.array([start]), grid)))
    return grid

def truncation(grid, inf, sup):
    if inf == None:
        inf = grid[0]
    if sup == None:
        sup = grid[-1]
    grid = grid[np.logical_and(grid>=inf, grid<=sup)]
    return grid

def stretch(grid, start, end):
    grid = start + (grid-grid[0])/(grid[-1]-grid[0])*(end-start)
    return grid

@jit(nopython=True, nogil=True, cache=True, parallel=True)
# @guvectorize([(float64, float64[:], float64[:], float64[:], float64[:], float64[:],
#               float64[:], float64[:], float64[:], float64[:], float64[:, :], float64[:],
#               float64[:])], target='parallel')
def corrected_jacobian_numba_simple(age_top_sigma,
                             accu,
                             depth, depth_inter, depth_mid,
                             age, age_model,
                             agedens, icelayerthick,
                             corr_a,
                             chol_a,
                             sigmap_corr_a,
                             corr_a_age):

    age_jac = np.zeros((1+len(corr_a), len(age)))
    age_jac[0, :] = age_top_sigma * np.ones(len(age))       
        
    for i in prange(len(corr_a)):

        corr_a_vec = np.zeros(len(corr_a))
        corr_a_vec[i] = 1.
    #Accu
        corr_vec = np.dot(chol_a, corr_a_vec)*sigmap_corr_a
        toto = np.interp((age_model[:-1]+age_model[1:])/2,
                                  corr_a_age, corr_vec)
        agedens_vec = - toto * agedens

    #Ice age
        age_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*agedens_vec)))
        age_jac[1+i,:] = age_vec
            

    return age_jac


@njit(nopython=True, nogil=True, cache=True, parallel=True)
def corrected_jacobian_numba_simple_full(age_top_sigma,
                              accu,
                              depth, depth_inter, depth_mid,
                              age, age_model,
                              agedens, icelayerthick,
                              corr_a,
                              chol_a, 
                              sigmap_corr_a, 
                              corr_a_age):

    accu_jac = np.zeros((1+len(corr_a), len(accu)))
    age_jac = np.zeros((1+len(corr_a), len(age)))
    age_jac[1, :] = age_top_sigma * np.ones(len(age))
        
    for i in prange(len(corr_a)):

        corr_a_vec = np.zeros(len(corr_a))
        corr_a_vec[i] = 1.
    #Accu
        corr_vec = np.dot(chol_a, corr_a_vec)*sigmap_corr_a
        toto = np.interp((age_model[:-1]+age_model[1:])/2,
                                  corr_a_age, corr_vec)
        agedens_vec = - toto * agedens
        accu_vec =  toto * accu
        accu_jac[1+i, :] = accu_vec

    #Ice age
        age_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*agedens_vec)))
        age_jac[1+i, :] = age_vec


    return accu_jac, age_jac


@njit(nopython=True, nogil=True, cache=True)
def corrected_jacobian_numba_icecore(age_top_sigma,
                             accu, tau, lid, dens, dens_firn,
                             depth, depth_inter, depth_mid,
                             age, airage, age_model, airage_model,
                             ice_equiv_depth,
                             agedens, icelayerthick,
                             corr_a, corr_tau, corr_lid,
                             chol_a, chol_tau, chol_lid,
                             sigmap_corr_a, sigmap_corr_tau, sigmap_corr_lid,
                             corr_a_age, corr_tau_depth, corr_lid_age):

    age_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(age)))
    age_jac[0, :] = age_top_sigma * np.ones(len(age))       
    airage_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(airage)))
    airage_jac[0, :] = age_top_sigma * np.ones(len(airage))       
    delta_depth_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(depth)))
#    delta_depth[0, :] = age_top_sigma * np.zeros(len(depth))  # useless

    for i in prange(len(corr_a)):

        corr_a_vec = np.zeros(len(corr_a))
        corr_a_vec[i] = 1.
    #Accu
        corr_vec = np.dot(chol_a, corr_a_vec)*sigmap_corr_a
        toto = np.interp((age_model[:-1]+age_model[1:])/2,
                                  corr_a_age, corr_vec)
        agedens_vec = - toto * agedens

    #Ice age
        age_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*agedens_vec)))
        age_jac[1+i,:] = age_vec

    #Air age
        airage_vec = np.interp(ice_equiv_depth, depth, age_vec)
        airage_jac[1+i, :] = airage_vec
        # delta_depth_vec = np.zeros_like(depth)
        # delta_depth_jac[1+i, :] = delta_depth_vec
            

        
    for i in prange(len(corr_tau)):
                        
        corr_tau_vec = np.zeros(len(corr_tau))
        corr_tau_vec[i] = 1.
        corr_vec = np.dot(chol_tau, corr_tau_vec)*sigmap_corr_tau
        tata = np.interp(depth_mid, corr_tau_depth, corr_vec)
        agedens_vec = -tata * agedens                
        age_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*agedens_vec)))
        age_jac[1+len(corr_a)+i, :] = age_vec
        
        thin_vec = -tata * dens/tau
        udepth_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*thin_vec)))
        delta_depth_vec = - np.interp(ice_equiv_depth, depth_mid,
                                      tau/dens) * (udepth_vec - \
                                    np.interp(ice_equiv_depth, depth, udepth_vec))
        airage_vec = np.interp(ice_equiv_depth, depth, age_vec) \
                        - np.interp(ice_equiv_depth, depth_mid, agedens) * \
                        delta_depth_vec
        airage_jac[1+len(corr_a)+i, :] = airage_vec
        delta_depth_jac[1+len(corr_a)+i, :] = delta_depth_vec

        #To be continued...                

    for i in prange(len(corr_lid)):

        age_vec = np.zeros_like(depth)
        age_jac[1+len(corr_a)+len(corr_tau)+i, :] = age_vec

        corr_lid_vec = np.zeros(len(corr_lid))
        corr_lid_vec[i] = 1.
        corr_vec = np.dot(chol_lid, corr_lid_vec)*sigmap_corr_lid
        lid_vec = np.interp(airage_model, corr_lid_age, corr_vec) * lid
        delta_depth_vec = dens_firn * lid_vec * \
                            np.interp(ice_equiv_depth, depth_mid, 
                                      tau/dens)
        airage_vec = - np.interp(ice_equiv_depth, depth_mid, 
                                  agedens) * delta_depth_vec
        airage_jac[1+len(corr_a)+len(corr_tau)+i, :] = airage_vec
        delta_depth_jac[1+len(corr_a)+len(corr_tau)+i, :] = delta_depth_vec

    return airage_jac, delta_depth_jac, age_jac


@njit(nopython=True, nogil=True, cache=True)
def corrected_jacobian_numba_icecore_full(age_top_sigma,
                              accu, tau, lid, dens, dens_firn,
                              depth, depth_inter, depth_mid,
                              age, airage, age_model, airage_model,
                              ice_equiv_depth,
                              agedens, icelayerthick,
                              corr_a, corr_tau, corr_lid,
                              chol_a, chol_tau, chol_lid, 
                              sigmap_corr_a, sigmap_corr_tau, sigmap_corr_lid, 
                              corr_a_age, corr_tau_depth, corr_lid_age):

    accu_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(accu)))
    age_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(age)))
    age_jac[1, :] = age_top_sigma * np.ones(len(age))
    airage_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(airage)))
    airage_jac[1, :] = age_top_sigma * np.ones(len(airage))
    delta_depth_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(depth)))
#        delta_depth_jac[1, :] = np.zeros(len(depth))
    icelayerthick_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(icelayerthick)))
#        icelayerthick_jac[1, :] = np.zeros(len(icelayerthick))
    tau_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(tau)))
    lid_jac = np.zeros((1+len(corr_a)+len(corr_tau)+len(corr_lid), len(lid)))
        
    for i in prange(len(corr_a)):

        corr_a_vec = np.zeros(len(corr_a))
        corr_a_vec[i] = 1.
    #Accu
        corr_vec = np.dot(chol_a, corr_a_vec)*sigmap_corr_a
        toto = np.interp((age_model[:-1]+age_model[1:])/2,
                                  corr_a_age, corr_vec)
        agedens_vec = - toto * agedens
        accu_vec =  toto * accu
        accu_jac[1+i, :] = accu_vec

    #Ice age
        age_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*agedens_vec)))
        age_jac[1+i, :] = age_vec

    #Air age
        airage_vec = np.interp(ice_equiv_depth, depth, age_vec)
        airage_jac[1+i, :] = airage_vec
        # delta_depth_vec = np.zeros_like(depth)
        # delta_depth_jac[1+i, :] = delta_depth_vec
        icelayerthick_vec = toto * icelayerthick
        icelayerthick_jac[1+i, :] = icelayerthick_vec
        # tau_vec = np.zeros_like(tau)
        # tau_jac[1+i, :] = tau_vec
        # lid_vec = np.zeros_like(lid)
        # lid_jac[1+i, :] = lid_vec
            

    
    for i in prange(len(corr_tau)):
                        
        corr_tau_vec = np.zeros(len(corr_tau))
        corr_tau_vec[i] = 1.
        corr_vec = np.dot(chol_tau, corr_tau_vec)*sigmap_corr_tau
        tata = np.interp(depth_mid, corr_tau_depth, corr_vec)
        agedens_vec = -tata * agedens                
        age_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*agedens_vec)))
        age_jac[1+len(corr_a)+i, :] = age_vec
        
        thin_vec = -tata * dens/tau
        udepth_vec = np.cumsum(np.concatenate((np.array([0]), depth_inter*thin_vec)))
        delta_depth_vec = - np.interp(ice_equiv_depth, depth_mid,
                                      tau/dens) * (udepth_vec - \
                                    np.interp(ice_equiv_depth, depth, udepth_vec))
        airage_vec = np.interp(ice_equiv_depth, depth, age_vec) \
                        - np.interp(ice_equiv_depth, depth_mid, agedens) * \
                        delta_depth_vec
        airage_jac[1+len(corr_a)+i, :] = airage_vec
        delta_depth_jac[1+len(corr_a)+i, :] = delta_depth_vec
        # accu_vec = np.zeros_like(accu)
        # accu_jac[1+len(corr_a)+i, :] = accu_vec
        icelayerthick_vec = tata * icelayerthick
        icelayerthick_jac[1+len(corr_a)+i, :] = icelayerthick_vec
        tau_vec = tata * tau
        tau_jac[1+len(corr_a)+i, :] = tau_vec
        # lid_vec = np.zeros_like(lid)
        # lid_jac[1+len(corr_a)+i, :] = lid_vec

        #To be continued...                

    for i in prange(len(corr_lid)):
                        
        # age_vec = np.zeros_like(depth)
        # age_jac[1+len(corr_a)+len(corr_tau)+i, :] = age_vec

        corr_lid_vec = np.zeros(len(corr_lid))
        corr_lid_vec[i] = 1.
        corr_vec = np.dot(chol_lid, corr_lid_vec)*sigmap_corr_lid
        lid_vec = np.interp(airage_model, corr_lid_age, corr_vec) * lid
        delta_depth_vec = dens_firn * lid_vec * \
                            np.interp(ice_equiv_depth, depth_mid, 
                                      tau/dens)
        airage_vec = - np.interp(ice_equiv_depth, depth_mid, 
                                  agedens) * delta_depth_vec
        airage_jac[1+len(corr_a)+len(corr_tau)+i, :] = airage_vec
        delta_depth_jac[1+len(corr_a)+len(corr_tau)+i, :] = delta_depth_vec
        # accu_vec = np.zeros_like(accu)
        # accu_jac[1+len(corr_a)+len(corr_tau)+i, :] = accu_vec
        # icelayerthick_vec = np.zeros_like(icelayerthick)
        # icelayerthick_jac[1+len(corr_a)+len(corr_tau)+i, :] = icelayerthick_vec
        tau_vec = np.zeros_like(tau)
        tau_jac[1+len(corr_a)+len(corr_tau)+i, :] = tau_vec
        lid_jac[1+len(corr_a)+len(corr_tau)+i, :] = lid_vec

    return accu_jac, airage_jac, delta_depth_jac, icelayerthick_jac,\
        tau_jac, lid_jac, age_jac
