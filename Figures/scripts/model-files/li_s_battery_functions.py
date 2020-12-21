# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 08:26:06 2019

@author: dkorff
"""

import cantera as ct
import numpy as np
from math import pi, tanh
from li_s_battery_inputs import inputs
from li_s_battery_init import cathode as cat

def dst(s1, s2, D_eff, dy1, dy2):
    F = ct.faraday; R = ct.gas_constant; T = inputs.T
    
    dyInv = 1/(0.5*(dy1 + dy2))
    w1 = dy2/(dy1 + dy2); w2 = dy1/(dy1 + dy2)
    C_k = (w1*s1['C_k'] + w2*s2['C_k'])  #*0.5
    z_k = inputs.z_k_el
    
    N_io = np.zeros_like(s1['C_k'])
    N_io = (-D_eff*(s2['C_k'] - s1['C_k'])*dyInv
            -D_eff*C_k*(z_k*F/R/T)*(s2['phi_el'] - s1['phi_el'])*dyInv)
    
#    N_io = N_io*[1, 1, 1, 0, 0, 0, 0, 0, 0]
    i_io = np.dot(N_io, z_k)*F
    
    return N_io, i_io

"""========================================================================="""

def set_state(SV, offset, ptr):
    
    state = {}
    
    state['phi_ed'] = SV[offset + ptr['phi_ed']]
    state['phi_dl'] = SV[offset + ptr['phi_dl']]
    state['phi_el'] = SV[offset + ptr['phi_ed']] - SV[offset + ptr['phi_dl']]
    state['C_tot'] = sum(SV[offset + ptr['rho_k_el']])
    state['C_k'] = SV[offset + ptr['rho_k_el']]
    state['X_k'] = SV[offset + ptr['rho_k_el']]/sum(SV[offset + ptr['rho_k_el']])
    
#    state['np_S'] = SV[offset + ptr['np_S8']]
#    state['np_L'] = SV[offset + ptr['np_Li2S']]
#    state['eps_S8'] = max(SV[offset + ptr['eps_S8']], 1e-25)
#    state['eps_Li2S'] = max(SV[offset + ptr['eps_Li2S']], 1e-25)
#    
#    state['A_S'] = 3*state['eps_S8']**(2/3)/(3*cat.V_0/2/pi/state['np_S'])**(1/3)
#    state['A_L'] = 3*state['eps_Li2S']/(3*state['eps_Li2S']*cat.V_0/2/pi/state['np_L'])**(1/3)
#    
#    state['r_S'] = 3*state['eps_S8']/state['A_S']
#    state['r_L'] = 3*state['eps_Li2S']/state['A_L']
#    
#    state['A_C'] = inputs.A_C_0 \
#                 - (pi*state['np_S']*state['r_S']**2)/cat.V_0 \
#                 - (pi*state['np_L']*state['r_L']**2)/cat.V_0
    
    return state

"""========================================================================="""

def set_geom(SV, offset, ptr):
    geom = {}
    
    geom['np_S'] = SV[offset + ptr['np_S8']]
    geom['np_L'] = SV[offset + ptr['np_Li2S']]
    geom['eps_S'] = max(SV[offset + ptr['eps_S8']], cat.eps_cutoff)
    geom['eps_L'] = max(SV[offset + ptr['eps_Li2S']], cat.eps_cutoff)
    geom['eps_el'] = 1 - cat.eps_C_0 - geom['eps_S'] - geom['eps_L']
    
    geom['A_S'] = 2*pi*geom['np_S']*(3*geom['eps_S']/2/geom['np_S']/pi)**(2/3)
    geom['A_L'] = 2*pi*geom['np_L']*(3*geom['eps_L']/2/geom['np_L']/pi)**(2/3)
    
    geom['r_S'] = 3*geom['eps_S']/geom['A_S']
    geom['r_L'] = 3*geom['eps_L']/geom['A_L']
    geom['L_tpb'] = 3*geom['eps_L']/(geom['r_L']**2)
    geom['A_C'] = cat.A_C_0 - (pi*geom['np_S']*geom['r_S']**2) \
                            - (pi*geom['np_L']*geom['r_L']**2)
    
    return geom

def set_rxn(geom, C_el_s, S_el_s, L_el_s, Li2S_tpb, sulfur, elyte, Li2S, conductor):
    sdot = {}
    
    sdot_C = C_el_s.get_net_production_rates(elyte)
    R_C = sdot_C*geom['A_C']
    
    mult = tanh(geom['eps_S']/cat.eps_dropoff)
    sdot['S8'] = S_el_s.get_creation_rates(sulfur) - \
            mult*S_el_s.get_destruction_rates(sulfur)
            
    sdot_S = S_el_s.get_net_production_rates(elyte)
    R_S = sdot_S*geom['A_S']
    
    mult = tanh(geom['eps_L']/cat.eps_dropoff)
    sdot['Li2S'] = L_el_s.get_creation_rates(Li2S) - \
              mult*L_el_s.get_destruction_rates(Li2S)
    sdot['tpb'] = Li2S_tpb.get_creation_rates(Li2S) - \
             mult*Li2S_tpb.get_destruction_rates(Li2S)
              
    sdot_L = L_el_s.get_net_production_rates(elyte)
    sdot_tpb_el = mult*Li2S_tpb.get_creation_rates(elyte) - \
                - Li2S_tpb.get_destruction_rates(elyte)
    R_L = sdot_L*geom['A_L'] + sdot_tpb_el*geom['L_tpb']
    R_net = R_C + R_S + R_L
    
    sdot['i_C'] = C_el_s.get_net_production_rates(conductor)*geom['A_C'] + \
                  Li2S_tpb.get_net_production_rates(conductor)*geom['L_tpb']
                  
    return sdot, R_net

"""========================================================================="""

def set_state_sep(SV, offset, ptr):
    
    state = {}
    
    state['phi_el'] = SV[offset + ptr['phi']]
    state['C_tot'] = sum(SV[offset + ptr['rho_k_el']])
    state['C_k'] = SV[offset + ptr['rho_k_el']]
    state['X_k'] = SV[offset + ptr['rho_k_el']]/sum(SV[offset + ptr['rho_k_el']])
    
    return state