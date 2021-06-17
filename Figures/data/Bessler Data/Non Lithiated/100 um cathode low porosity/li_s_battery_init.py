# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:39:12 2019

@author: dkorff

This is the object initialization file for the Li-s model. It imports input
values from li_s_battery_inputs.py and initializes the necessary objects to
run the simulations
"""

import numpy as np
import cantera as ct
import importlib
from math import pi

from li_s_battery_inputs import inputs

"Import cantera objects - this step is the same regardless of test type"
elyte_obj = ct.Solution(inputs.ctifile, inputs.elyte_phase)
sulfur_obj = ct.Solution(inputs.ctifile, inputs.cat_phase1)
Li2S_obj = ct.Solution(inputs.ctifile, inputs.cat_phase2)
carbon_obj = ct.Solution(inputs.ctifile, inputs.cat_phase3)
conductor_obj = ct.Solution(inputs.ctifile, inputs.metal_phase)
lithium_obj = ct.Solution(inputs.ctifile, inputs.an_phase)

sulfur_el_s = ct.Interface(inputs.ctifile, inputs.sulfur_elyte_phase,
                             [sulfur_obj, elyte_obj, conductor_obj])
Li2S_el_s = ct.Interface(inputs.ctifile, inputs.Li2S_elyte_phase,
                             [Li2S_obj, elyte_obj, conductor_obj])
carbon_el_s = ct.Interface(inputs.ctifile, inputs.graphite_elyte_phase,
                             [carbon_obj, elyte_obj, conductor_obj])
lithium_el_s = ct.Interface(inputs.ctifile, inputs.anode_elyte_phase,
                             [lithium_obj, elyte_obj, conductor_obj])
Li2S_tpb = ct.Interface(inputs.ctifile, 'tpb', [elyte_obj, Li2S_obj, conductor_obj])

elyte_obj.electric_potential = inputs.Phi_el_init
carbon_obj.electric_potential = inputs.Cell_voltage
conductor_obj.electric_potential = inputs.Cell_voltage
#print('k_f =', lithium_el_s.forward_rate_constants)
#print('k_r =', lithium_el_s.reverse_rate_constants)
#print('dG =', Li2S_el_s.delta_standard_gibbs)
#print('E0 =', -lithium_el_s.delta_standard_gibbs/ct.faraday)

#print(sulfur_el_s.forward_rate_constants, '\n', sulfur_el_s.reverse_rate_constants)
#print(carbon_el_s.forward_rate_constants, '\n', carbon_el_s.reverse_rate_constants)
#print(Li2S_el_s.forward_rate_constants, '\n', Li2S_el_s.reverse_rate_constants)

if 'cascade' or 'Bessler' in inputs.ctifile:
    dG0_an = lithium_el_s.delta_standard_gibbs
    E0_an = -dG0_an/ct.faraday
    print(E0_an)
    dG0 = carbon_el_s.delta_standard_gibbs
    E0 = dG0/ct.faraday
    print(-E0)
#    dG0_Li2S = Li2S_el_s.delta_standard_gibbs
#    print(dG0_Li2S)

if hasattr(inputs, 'C_k_el_0'):
    elyte_obj.X = inputs.C_k_el_0/np.sum(inputs.C_k_el_0)

bc_class = getattr(inputs, 'test_type')

F = ct.faraday       
       
"============================================================================="

class cathode():
    # Set a flag to let the solver know whether to implement this class
    flag = inputs.flag_cathode
    F = ct.faraday
    
    # Number of nodes in the y-direction
    npoints = inputs.npoints_cathode
        
    # Number of state variables per node
    nVars = 2 + elyte_obj.n_species + 4
    
    # Pointers
    ptr = {}
    ptr['iFar'] = elyte_obj.species_index(inputs.Li_species_elyte)
    ptr['eps_S8'] = 0
    ptr['eps_Li2S'] = 1
    ptr['rho_k_el'] = 2 + np.arange(0, elyte_obj.n_species)
    ptr['phi_dl'] = ptr['rho_k_el'][-1] + 1
    ptr['phi_ed'] = ptr['rho_k_el'][-1] + 2
    ptr['np_S8'] = ptr['rho_k_el'][-1] + 3
    ptr['np_Li2S'] = ptr['rho_k_el'][-1] + 4
    
    nSV = npoints*nVars
    offsets = np.arange(0, int(nSV), int(nVars))
    
    ptr_vec = {}
    ptr_vec['eps_S8']   = ptr['eps_S8']   + offsets
    ptr_vec['eps_Li2S'] = ptr['eps_Li2S'] + offsets
    ptr_vec['rho_k_el'] = ptr['rho_k_el']
    for i in offsets[1:]:
        ptr_vec['rho_k_el'] = np.hstack((ptr_vec['rho_k_el'],i+ptr['rho_k_el']))
    ptr_vec['phi_dl']  = ptr['phi_dl'] + offsets
    ptr_vec['phi_ed']  = ptr['phi_ed']   + offsets
    ptr_vec['np_S8']   = ptr['np_S8']   + offsets
    ptr_vec['np_Li2S'] = ptr['np_Li2S'] + offsets
    
    # Store parameters as class attributes
    T = inputs.T
    C_dl = inputs.C_dl_cat
    
    # Geometric parameters
    tau = inputs.tau_cat
    r_p = inputs.r_p_cat
    d_p = inputs.d_p_cat
    dyInv = npoints/inputs.H_cat
    dy = inputs.H_cat/npoints
    H = inputs.H_cat
    V_0 = inputs.H_cat*inputs.A_cat
    
    
    if inputs.sulfur_method == 'bulk':
        m_S = inputs.m_S_0/inputs.A_cat
        m_S_0 = inputs.m_S_0
    elif inputs.sulfur_method == 'loading':
        m_S = inputs.m_S_0
        m_S_0 = inputs.m_S_0*inputs.A_cat
        
        
    omega_S = inputs.pct_w_S8_0/inputs.A_cat
    omega_C = inputs.pct_w_C_0/inputs.A_cat
    rho_S = sulfur_obj.density_mass
    rho_C = carbon_obj.density_mass
    m_solid = m_S/omega_S
    
#    eps_S_0 = 0.16  
    eps_S_0 = m_S/rho_S/H
#    eps_C_0 = 0.062 
    eps_C_0 = 5*0.056  #m_solid*omega_C/rho_C/H
    print('Eps_C_0 =', eps_C_0)
#    eps_L_0 = 1e-4; 
    eps_L_0 = 1e-5
    
#    if inputs.mech == 'Bessler-Dennis':
#        m_S_0 = eps_S_0*H*inputs.A_cat*sulfur_obj.density_mass
    
#    A_S_0 = 1e5  
    A_S_0 = 2*pi*inputs.np_S8_init*(3*eps_S_0/2/inputs.np_S8_init/pi)**(2/3)
#    A_L_0 = 1e5  
    A_L_0 = 2*pi*inputs.np_Li2S_init*(3*eps_L_0/2/inputs.np_Li2S_init/pi)**(2/3)
    
    r_S_0 = 3*eps_S_0/A_S_0
    r_L_0 = 3*eps_L_0/A_L_0
#    print(r_S_0, r_L_0)

    A_C_0 = inputs.A_C_0
    print('A_S =', A_S_0)
    print('A_L =', A_L_0)
    print('A_C =', A_C_0 - (pi*inputs.np_S8_init*r_S_0**2) - (pi*inputs.np_Li2S_init*r_L_0**2))
    
    eps_el_0 = 1 - eps_S_0 - eps_C_0 - eps_L_0
    eps_pore = 1 - eps_C_0
    
    
    print('Porosity =', eps_el_0)
    
    m_el = H*eps_el_0*elyte_obj.density_mass
    m_sulfur = H*eps_S_0*sulfur_obj.density_mass
    m_carbon = H*eps_C_0*carbon_obj.density_mass
    m_L = H*eps_L_0*Li2S_obj.density_mass
    m_cat = inputs.A_cat*(m_el + m_sulfur + m_carbon + m_L)
    
    n_S_atoms = np.zeros([len(elyte_obj.species_names)])
    for i, species in enumerate(elyte_obj.species_names):
        if elyte_obj.n_atoms(species, 'S') and i != 2:
            n_S_atoms[i] = elyte_obj.n_atoms(species, 'S')
        else:
            n_S_atoms[i] = 0
            
    S_atoms_bool = np.zeros_like(n_S_atoms)
    for i in np.arange(len(n_S_atoms)):
        if n_S_atoms[i] == 0:
            S_atoms_bool[i] = 0
        elif n_S_atoms[i] != 0 and i != 2:
            S_atoms_bool[i] = 1
            
    n_S_species = np.count_nonzero(S_atoms_bool)
        
    n_S_0 = eps_el_0*H*np.dot(n_S_atoms, inputs.C_k_el_0) \
          + 8*sulfur_obj.density_mole*eps_S_0*H \
          + Li2S_obj.density_mole*eps_L_0*H
              
#    W_S_k = elyte_obj.molecular_weights*S_atoms_bool # Old method
    W_S = sulfur_obj.molecular_weights/sulfur_obj.n_atoms(sulfur_obj.species_names[0], 'S')
    cap_weights = np.array([1, 7/8, 0.8333, 0.75, 0.5, 0])
    W_S_k = elyte_obj.molecular_weights[3:]
    
    m_S_el = inputs.A_cat*eps_el_0*H*W_S*np.dot(n_S_atoms, inputs.C_k_el_0)
    m_S_el_an = inputs.A_cat*(1 - inputs.epsilon_an)*inputs.H_an*W_S*np.dot(n_S_atoms, inputs.C_k_el_0)
    m_S_el_sep = inputs.A_cat*(1 - inputs.epsilon_sep)*inputs.H_elyte*W_S*np.dot(n_S_atoms, inputs.C_k_el_0)
    m_S_tot_0 = m_S_0 + m_S_el + m_S_el_an + m_S_el_sep 
    
    V_elyte = inputs.A_cat*(inputs.H_an*(1-inputs.epsilon_an) +
                            inputs.H_elyte*(1-inputs.epsilon_sep) +
                            inputs.H_cat*eps_el_0)
    print('Elyte/sulfur ratio ', 1e3*V_elyte/m_S_tot_0)
    
    x = np.copy(n_S_atoms)
    x[5:] = (x[5:] - 1)
    oneC = 1675*(m_S_tot_0)/inputs.A_cat
    print('solid sulfur =', m_S_0/inputs.A_cat)
    
    def get_i_ext():
        return cathode.i_ext
    
    def set_i_ext(value):
        cathode.i_ext = value
        
    nucleation_flag = np.zeros((inputs.npoints_cathode, 1))
    np_L = inputs.np_Li2S_init*np.ones((inputs.npoints_cathode, 1))
            
    # Calculate the actual current density. 
    i_ext_amp = -inputs.C_rate*oneC
    print('External current =', i_ext_amp)
    
    sigma_eff = inputs.sigma_cat*eps_C_0/tau**3
        
    D_el = inputs.D_Li_el
    bruggeman = 1.5
    
    eps_cutoff = 1e-15
    eps_S8_cutoff = 1e-5
    eps_dropoff = 1e-10
    
#    z_k_el = elyte_obj.species_charges
    
    def get_tflag():
        return cathode.t_flag
    
    def set_tflag(value):
        cathode.t_flag = value
        
    def set_tags(value):
        cathode.tags = value
        
    A_C_vec = np.array([])
    nucl_thresh = 1e-2
    
"============================================================================="        
        
class sep():
    # Set a flag to let the solver know whether to implement this class
    flag = inputs.flag_sep
    
    # Number of nodes in the y-direction
    npoints = inputs.npoints_sep
    
    # Number of variables per node
    nVars = 1 + elyte_obj.n_species
    
    H = inputs.H_elyte  # Separator thickness [m]
    
    tau = inputs.tau_sep  # Tortuosity of separator
    
    # Geometric parameters
    epsilon = inputs.epsilon_sep  # Volume fraction of separator material [-]
    epsilon_el = 1 - epsilon      # Volume fraction of electrolyte [-]
    dyInv = npoints/H             # Inverse of y-direction discretization [1/m]
    dy = H/npoints
    sep_density_mass = 940        # HDPE density [kg/m^3]
    # Mobility of electrolyte species
    u_Li_el = inputs.D_Li_el*epsilon_el/ct.gas_constant/inputs.T/tau**3
    
    ptr = {}
    ptr['rho_k_el'] = np.arange(0, elyte_obj.n_species)
    ptr['phi'] = elyte_obj.n_species
    
    ptr_vec = {}
    ptr_vec['rho_k_el'] = cathode.nSV + ptr['rho_k_el']
    ptr_vec['phi'] = cathode.nSV + ptr['phi']
    
    for i in np.arange(1, npoints):
        ptr_vec['rho_k_el'] = np.append(ptr_vec['rho_k_el'], 
                                      cathode.nSV + ptr['rho_k_el'] + i*nVars)
        ptr_vec['phi'] = np.append(ptr_vec['phi'], 
                                   cathode.nSV + ptr['phi'] + i*nVars)
        
    # Set the length of the solution vector for the separator
    nSV = npoints*nVars
    
    D_el = inputs.D_Li_el*epsilon_el**(cathode.bruggeman)
    
    offsets = np.arange(int(cathode.nSV), int(cathode.nSV) + int(nSV), int(nVars))
    
    n_S_0 = epsilon_el*H*np.dot(cathode.n_S_atoms, inputs.C_k_el_0)
    
    m_HDPE = H*epsilon*sep_density_mass
    m_el = H*epsilon_el*elyte_obj.density_mass
    m_sep = inputs.A_cat*(m_HDPE + m_el)
    
#    z_k_el = elyte_obj.species_charges
    
"============================================================================="

class anode():
    flag = inputs.flag_anode
    
    npoints = inputs.npoints_anode
    
    nVars = 2 + elyte_obj.n_species
    
    # Pointers
    ptr = {}
    ptr['iFar'] = elyte_obj.species_index(inputs.Li_species_elyte)
    
    ptr['rho_k_el'] = np.arange(0, elyte_obj.n_species)
    ptr['phi_dl'] = ptr['rho_k_el'][-1] + 1
    ptr['phi_ed'] = ptr['rho_k_el'][-1] + 2
    
    ptr_vec = {}
    ptr_vec['rho_k_el'] = cathode.nSV + sep.nSV + ptr['rho_k_el']
    
    for i in np.arange(1, npoints):
        ptr_vec['rho_k_el'] = np.append(ptr_vec['rho_k_el'],
                                       cathode.nSV + sep.nSV + ptr['rho_k_el'] + i*nVars)
    
    # Set length of solution vector for anode
    nSV = npoints*nVars
    offsets = np.arange(int(cathode.nSV + sep.nSV), 
                        int(cathode.nSV + sep.nSV) + int(nSV), int(nVars))
    
    # Geometric parameters
    eps_el = 1 - inputs.epsilon_an
    tau = inputs.tau_an
    r_p = inputs.r_p_cat
    dyInv = npoints/inputs.H_an
    dy = inputs.H_an/npoints
    H = inputs.H_an
    
    H_el = H*eps_el
    dy_el = H_el/1
    dyInv_el = 1/H_el
    
    C_dl = inputs.C_dl_an
    A_Li = 1.5  #1/H
    sigma_eff = inputs.sigma_an*inputs.epsilon_an/tau**3
    
    u_Li_el = inputs.D_Li_el*eps_el/tau**3
    
    D_el = inputs.D_Li_el*eps_el**(1.)/tau**3
    
    n_S_0 = eps_el*H*np.dot(cathode.n_S_atoms, inputs.C_k_el_0)
    
    m_Li = H*inputs.epsilon_an*lithium_obj.density_mass
    m_el = H*eps_el*elyte_obj.density_mass
    m_an = inputs.A_cat*(m_Li + m_el)
    m_bat = cathode.m_cat + sep.m_sep + m_an
    
#    z_k_el = elyte_obj.species_charges
        
"============================================================================="

class sol_init():
    
    # Initialize solution vector 
    SV_0 = np.zeros([anode.nSV + sep.nSV + cathode.nSV])
    
    # Set up algebraic variable vector
    algvar = np.zeros_like(SV_0)
     
    # Cathode
    offsets = cathode.offsets
    ptr = cathode.ptr
    for j in np.arange(0, cathode.npoints):
        SV_0[offsets[j] + ptr['eps_S8']] = cathode.eps_S_0
        algvar[offsets[j] + ptr['eps_S8']] = 1
        
        SV_0[offsets[j] + ptr['eps_Li2S']] = cathode.eps_L_0
        algvar[offsets[j] + ptr['eps_Li2S']] = 1
        
        SV_0[offsets[j] + ptr['rho_k_el']] = inputs.C_k_el_0
        algvar[offsets[j] + ptr['rho_k_el']] = 1
        
        SV_0[offsets[j]+ptr['phi_dl']] = inputs.Cell_voltage - inputs.Phi_el_init
        algvar[offsets[j] + ptr['phi_dl']] = 1
                                           
        SV_0[offsets[j]+ptr['phi_ed']] = inputs.Cell_voltage
        
        SV_0[offsets[j]+ptr['np_S8']]=inputs.np_S8_init
#        algvar[offsets[j] + ptr['np_S8']] = 1
        
        SV_0[offsets[j]+ptr['np_Li2S']] = inputs.np_Li2S_init
#        algvar[offsets[j] + ptr['np_Li2S']] = 1
     
    # Separator
    offsets = sep.offsets
    ptr = sep.ptr
    for j in np.arange(0, sep.npoints):
        
        SV_0[offsets[j] + ptr['rho_k_el']] = inputs.C_k_el_0
        algvar[offsets[j] + ptr['rho_k_el']] = 1
        
        SV_0[offsets[j] + ptr['phi']] = inputs.Phi_el_init
        algvar[offsets[j] + ptr['phi']] = 1
     
    # Anode
    offsets = anode.offsets
    ptr = anode.ptr
    for j in np.arange(0, anode.npoints):
        SV_0[offsets[j] + ptr['rho_k_el']] = inputs.C_k_el_0
        algvar[offsets[j] + ptr['rho_k_el']] = 1
        
        SV_0[offsets[j] + ptr['phi_dl']] = inputs.Phi_an_init - inputs.Phi_el_init
        algvar[offsets[j] + ptr['phi_dl']] = 1
        
        SV_0[offsets[j] + ptr['phi_ed']] = inputs.Phi_an_init
    
                                           
"============================================================================="

print("Initialization check")