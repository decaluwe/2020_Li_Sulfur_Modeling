# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:02:23 2019

@author: dkorff

The code here and related modules models a lithium-sulfur battery (Li-S). It's
structured to allow a variety of test modes. Currently this includes constant
current dis/charge cycling and sinusoidal current impedence spectroscopy.

This module receives all inputs from the user and passes those inputs to an
initialization file that creates the necessary objects for the simulations
and then passes those objects to the solver file. The solver file will call
the appropriate residual function tied to the numerical solver. Post-processing
will be handled by a standalone module
"""

import numpy as np
from math import pi, exp

class inputs():
    
    # The following flags set which components to include in the simulation 
    flag_anode = 1
    flag_sep = 1
    flag_cathode = 1
    
    n_comps = flag_anode + flag_sep + flag_cathode
    
    # Set number of discretized nodes in each component's y-direction
    npoints_anode = 1*flag_anode
    npoints_sep = 5*flag_sep
    npoints_cathode = 20*flag_cathode
    
    # Set number of discretized shells in each particle    
    flag_req = 0
    n_cycles = 1
    
    """Plotting options"""
    # To turn on profile plotting set to 1
    flag_plot_profiles = 0
    
    # To plot voltage profile set to 1
    flag_potential = 1*flag_plot_profiles
    flag_electrode = 1*flag_plot_profiles
    flag_electrolyte = 1*flag_plot_profiles
    flag_capacity = 1*flag_plot_profiles
    
    # The C-rate is the rate of charge/discharge - how many charges/discharges
    #   can be carried out in 1 hour theoretically? This sets current density
    #   amplitude for impedence tests and external current for CC cycling
    C_rate = 1.0
#    C_rate = 1
    
    # Set the test type to run the model for. The following types are supported
    #   For constant external current dis/charge cycling test set to:
    #       test_type = cc_cycling
    #   For sinusoidal external current impedence spectroscopy set to:
    #       test_type = is_cycling
    test_type = 'cc_cycling'
    
    # Set temperature for isothermal testing
    T = 298.15  # [K]
    
    "Set up Cantera phase names and CTI file info"
#    ctifile = 'sulfur_cathode_cascade_Crate.cti'
#    ctifile = 'sulfur_cathode_cascade_lithiated.cti'
#    ctifile = 'Kuzmina.yml'
#    ctifile = 'Kuzmina3.yml'
#    ctifile = 'Assary.yml'
#    ctifile = 'Bessler_Dennis.yml'
#    ctifile = 'Bessler_Dennis_mod.yml'
    ctifile = 'Bessler_Dennis_lithiated.yml'
#    ctifile = 'Shriram.yml'
#    ctifile = 'Shriram_adjusted.yml'
    cat_phase1 = 'sulfur'
    cat_phase2 = 'lithium_sulfide'
    cat_phase3 = 'carbon'
    metal_phase = 'electron'
    elyte_phase = 'electrolyte'
    an_phase = 'lithium'
    
    sulfur_elyte_phase = 'sulfur_surf'
    graphite_elyte_phase = 'carbon_surf'
    Li2S_elyte_phase = 'lithium_sulfide_surf'
    tpb_phase = 'tpb'
    anode_elyte_phase = 'lithium_surf'
    
    Li_species_elyte = 'Li+(e)'
    Max_sulfide = 'S8(e)'
    
#    # Set initial SOC 
#    SOC_0 = 0.1
    
    # Set initial potential values for anode, elyte, and cell
    Phi_an_init = 0.0
    if 'Bessler' in ctifile:
        Phi_el_init = 0
    else:
        Phi_el_init = 1.3
    Cell_voltage = 2.4

    # Cutoff values for charging and discharging of electrodes:
    Li_an_min = 0.01; Li_an_max = 1 - Li_an_min
    Li_cat_min = 0.01; Li_cat_max = 1 - Li_cat_min
    
    # Cell geometry
    H_cat = 100e-6               # Cathode thickness [m]
    r_C = H_cat/npoints_cathode/2
#    A_C_0 = 1.32e5  # Initial volume specific area of carbon [1/m]
    A_C_0 = 5*2e4
    
    # There are two options for providing sulfur loading. Input the value in
    #   [kg_sulfur/m^2] pre-calculated or enter the mass of sulfur and cell
    #   area separately in [kg_sulfur] and [m^2] respectively. Enter 'loading'
    #   or 'bulk' in the string >sulfur_method below.
    sulfur_method = 'loading'
    A_cat = 1.327e-4            # Cathode planar area [m^2]
    m_S_0 = 3.57e-2 #1.9e-2     # Initial total mass of sulfur in cathode [kg_S8] 2.5e-2
                    #1.957e-2   # if 'bulk' method chosen. Sulfur loading in
                    #2.53e-2    # [kg_S8/m^2] if 'loading' method chosen.
                    #3.57e-2 
                    
    # Initial number of nucleation sites per volume for solid phases. Eventually will
    #   use a nucleation theory.
    if 'cascade' or 'Bessler' in ctifile:
#        n = 620999563729888.0  #1520999563729888.0  #275771605112807.8
        n = 6e15*exp(2.4221*C_rate)  #5e13*exp(2.4221*C_rate)  
        print("Density for cascade")
    else:
        n = 6e13*exp(1.8966*C_rate)  #8e12*exp(1.7953*C_rate)
        print("Density for Assary or Kuzmina", n)
    
#    n = 5e13
#    np_S8_init = npoints_cathode*1000/H_cat/A_cat
    np_S8_init = 4521477015825  #npoints_cathode*3000/H_cat/A_cat # 1000000 Initial number of sulfur nucleation sites [n/m^3]
    np_Li2S_init = n   # Initial number of Li2S nucleation sites [n/m^3]
    
    # Weight percent of sulfur in the cathode per cathode volume, this assumes 
    #   the complementary phase is only the carbon structure - i.e. 40 wt% 
    #   sulfur means 60 wt% carbon.
    pct_w_S8_0 = 0.6  # Initial weight percent of sulfur in cathode [kg_S8/kg]
    pct_w_C_0 = 1 - pct_w_S8_0   # Initial weight percent of carbon in cathode [kg_C/kg]
    C_counter_n = (1.024 - 1.821e-4*2 - 3.314e-4*2 - 2.046e-5*2 - 
                    5.348e-10*2 - 8.456e-10*2)
    if 'Kuzmina' in ctifile:
        C_k_el_0 = np.array([1.023e1, 
                             1.024, 
                             1.024, 
                             1.943e-4, 
                             1.821e-4, 
                             3.314e-6, 
                             2.046e-8,
                             2.046e-10,
                             5.348e-13])
        z_k_el = np.array([0., 1., -1., 0., 0., 0., 0., 0, 0.])
        # DK: Need to rebuild with updated version of cantera to use `species_charges`
        #   then parameter will be set in the _init.py file
        mech = 'Kuzmina'
        print('Using Kuzmina')
    elif 'Assary' in ctifile:
        C_k_el_0 = np.array([1.023e1, 
                             1.024, 
                             1.024, 
                             1.943e-4, 
                             1.821e-4, 
                             3.314e-6, 
                             2.046e-8,
                             2.046e-10,
                             5.348e-16])
        z_k_el = np.array([0., 1., -1., 0., 0., 0., 0., 0, 0.])
        mech = 'Assary'
        print('Using Assary')
    elif 'cascade' in ctifile:
        C_counter_n = (1.024 - 1.821e-4*2 - 3.314e-6*2 - 
                        2.046e-6*2 - 2.046e-6*2 - 5.348e-6*2)
        if 'lithiated' in ctifile:
            C_counter_0 = 1.024
        else:
            C_counter_0 = C_counter_n
            
        C_k_el_0 = np.array([1.023e1, 
                             1.024, 
                             C_counter_0, 
                             1.943e-2, 
                             1.821e-4, 
                             3.314e-6, 
                             2.046e-6,
                             2.046e-6,
                             5.348e-6])
        
        if 'lithiated' in ctifile:
            mech = 'Cascade_li'
            z_k_el = np.array([0., 1., -1., 0., 0., 0., 0., 0., 0.])
        else:
            mech = 'Cascade'
            z_k_el = np.array([0., 1., -1., 0., -2., -2., -2., -2., -2.])
        print('Using cascade')
    elif 'Bessler' or 'Shriram' in ctifile:
        C_counter_n = 1.024 - 1.821e-4*2 - 3.314e-4*2 - 2.046e-5*2 - 5.348e-10*2 - 8.456e-13*2
        if 'lithiated' in ctifile:
            C_counter_0 = 1.024
        else:
            C_counter_0 = C_counter_n
        C_k_el_0 = np.array([1.023e1, 
                             1.024, 
                             C_counter_0, 
                             1.943e-2, 
                             1.821e-4, 
                             3.314e-4, 
                             2.046e-5,
                             5.348e-10,
                             8.456e-13])
        if 'lithiated' in ctifile:
            mech = 'Bessler_Li'
            z_k_el = np.array([0., 1., -1., 0., 0., 0., 0., 0., 0.])
        else:
            mech = 'Bessler'
            z_k_el = np.array([0., 1., -1., 0., -2., -2., -2., -2., -2.])
    elif 'Shriram' in ctifile:
        C_counter_n = 1.0
        C_k_el_0 = np.array([1.023e1,
                             1.00104,
                             1.0,
                             1.9e-2,
                             1.78e-4,
                             3.24e-4,
                             2.0e-5,
                             5.229e-10,
                             8.267e-13])
        z_k_el = np.array([0., 1., -1., 0., -2., -2., -2., -2., -2.])
        print('Using Shriram')
        mech = 'Shriram'
    
    "Cathode geometry and transport"
    # Anode geometry
    epsilon_carbon = 0.062      # Volume fraction of carbon in cathode [-]
    tau_cat = 1.6               # Tortuosity for cathode [-]
    r_p_cat = 5e-6              # Average pore radius [m]
    d_p_cat = 5e-6              # Average particle diameter [m]
    overlap_cat = 0.4           # Percentage of carbon particle overlapping with other
                                #   carbon particles. Reduces total carbon/elyte
                                #   surface area [-]
                        
    # Transport
    C_dl_cat = 1.5e-2    # Double-layer capacitance [F/m^2]
    sigma_cat = 60.0     # Bulk cathode electrical conductivity [S/m]
    
    "Anode geometry and transport"
    # Anode geometry
    epsilon_an = 0.95  #0.63    # Volume fraction of anode phase [-]
    tau_an = 1.6        # Tortuosity, assume equal values for carbon and elyte [-]
    r_p_an = 5e-6       # Average pore radius [m]
    d_p_an = 5e-6       # Average particle diameter for graphite [m]
    H_an = 100e-6        # Anode thickness [m]
    overlap_an = 0.4    # Percentage of anode particle overlapping with other
                        #   anode particles. Reduces total anode/elyte
                        #   surface area [-]
                        
    # Transport
    C_dl_an = 1.5e-1    # Double-layer capacitance [F/m^2]
    sigma_an = 75.0     # Bulk anode electrical conductivity [S/m]
    D_Li_an = 7.5e-16   # Bulk diffusion coefficient for Li in graphite [m^2/s]
    
    "Electrolyte/separator geometry and transport"
    H_elyte = 25e-6     # Separator thickness [m]
    
    # Elytespecies bulk diffusion coefficients and charges. 
    #   Assumes four component elyte: [TEGDME, Li+, TFSI-, S8(e), Li2S8, Li2S6
    #                                  Li2S4, Li2S3, Li2S2]
#    D_Li_el = np.array([1e-12, 1e-10, 4e-10, 1e-11, 6e-11, 6e-11, 1e-10,
#                        1e-10, 1e-10])
    D_Li_el = np.array([1e-12, 1e-10, 4e-10, 1e-11, 6e-11, 6e-11, 1e-10,
                        1e-10, 1e-10])
    
    epsilon_sep = 0.5   # Volume fraction of separator [-]
    tau_sep = 1.6       # Tortuosity of separator [-]
    sigma_sep = 50.0    # Bulk ionic conductivity of separator [S/m]
    
print("Inputs check")

if __name__ == "__main__":
    exec(open("li_s_battery_init.py").read())
    exec(open("li_s_battery_model.py").read())
    
    
    