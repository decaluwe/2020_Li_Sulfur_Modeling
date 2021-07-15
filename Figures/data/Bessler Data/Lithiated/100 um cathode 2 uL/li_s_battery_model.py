# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:30:28 2019

@author: dkorff

This module contains the set up and running of the simulation for the Li-S
model.
"""

import numpy as np
import pandas as pd
import time
import importlib
import cantera as ct
from matplotlib import pyplot as plt
import matplotlib
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem
from assimulo.exception import TerminateSimulation

from li_s_battery_inputs import inputs

from li_s_battery_init import anode as an
from li_s_battery_init import sep
from li_s_battery_init import cathode as cat
from li_s_battery_init import sol_init

from li_s_battery_post import label_columns
from li_s_battery_post import tag_strings
from li_s_battery_post import plot_sim
from li_s_battery_post import plot_meanPS  

from li_s_battery_init import elyte_obj as elyte

def main():
    
    res_class = eval(inputs.test_type)
    
#    plt.close('all')
    t_count = time.time()
        
    SV_0 = sol_init.SV_0
    SV_dot_0 = np.zeros_like(SV_0)
    t_0 = 0.
    t_f = 3600./inputs.C_rate  #63006.69900049 93633
    algvar = sol_init.algvar
    atol = np.ones_like(SV_0)*1e-3
    atol[cat.ptr_vec['eps_S8']] = 1e-22
    atol[cat.ptr_vec['eps_Li2S']] = 1e-15
    atol[cat.ptr_vec['rho_k_el']] = 1e-26 # 1e-19 for Bessler
#    atol[cat.ptr_vec['rho_k_el']][3::elyte.n_species] = 26
    rtol = 1e-4; sim_output = 50  # 3e-4 at 1C
     
    rtol_ch = 1e-6
    atol_ch = np.ones_like(SV_0)*1e-6
    
    atol_ch[cat.ptr_vec['eps_S8']] = 1e-15
    atol_ch[cat.ptr_vec['eps_Li2S']] = 1e-15
    atol_ch[cat.ptr_vec['rho_k_el']] = 1e-30
    
    rate_tag = str(inputs.C_rate)+"C"
#    if 'cascade' in inputs.ctifile:
#        ncols = 2 + inputs.flag_req
#    else:
    ncols = 1 + inputs.flag_req
    fig, axes = plt.subplots(sharey="row", figsize=(9,12), nrows=3, ncols = (ncols)*inputs.n_cycles, num=1)
    plt.subplots_adjust(wspace = 0.15, hspace = 0.4)
    fig.text(0.35, 0.85, rate_tag, fontsize=20, bbox=dict(facecolor='white', alpha = 0.5))
    
    # Set up user function to build figures based on inputs
    
    "----------Equilibration----------"
    
    print('\nEquilibrating...')

    # Set external current to 0 for equilibration
    cat.set_i_ext(0)
#    cat.nucleation_flag = 1
    
    # Create problem object
    bat_eq = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
    bat_eq.external_event_detection = True
    bat_eq.algvar = algvar
    
    # Create simulation object
    sim_eq = IDA(bat_eq)
    sim_eq.atol = atol
    sim_eq.rtol = rtol
    sim_eq.verbosity = sim_output
    sim_eq.make_consistent('IDA_YA_YDP_INIT')
    
    t_eq, SV_eq, SV_dot_eq = sim_eq.simulate(t_f)
    
    # Put solution into pandas dataframe with labeled columns
    SV_eq_df = label_columns(t_eq, SV_eq, an.npoints, sep.npoints, cat.npoints)
    
    # Obtain tag strings for dataframe columns
    tags = tag_strings(SV_eq_df)
#    plot_sim(tags, SV_eq_df, 'Equilibrating', 0, fig, axes)
    
    t_equilibrate = time.time() - t_count
    print("Equilibration time = ", t_equilibrate, '\n')
    
    print('Done equilibrating\n')
        
    for cycle_num in np.arange(0, inputs.n_cycles):
        "------------Discharging-------------"
    
        print('Discharging...')
#        cat.np_L = inputs.np_Li2S_init
        cat.nucleation_flag = np.zeros((inputs.npoints_cathode, 1))
        # New initial conditions from previous simulation
        if cycle_num == 0:
#            SV_0 = SV_0
            SV_0 = SV_eq[-1, :]
#            SV_dot_0 = SV_dot_0
            SV_dot_0 = SV_dot_eq[-1, :]
        else:   
            SV_0 = SV_0_cycle  
            SV_dot_0 = SV_dot_0_cycle  
        
        # Set external current
        cat.set_i_ext(cat.i_ext_amp)
        
        # Update problem instance initial conditions
        bat_dch = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
        bat_dch.external_event_detection = True
        bat_dch.algvar = algvar
            
        # Re-initialize simulation object
        sim_dch = IDA(bat_dch)
        sim_dch.atol = atol
        sim_dch.rtol = rtol
        sim_dch.maxh = 2
        sim_dch.inith = 1e-5  #1e-5 for bessler
        sim_dch.verbosity = sim_output
        sim_dch.make_consistent('IDA_YA_YDP_INIT')

        t_dch, SV_dch, SV_dot_dch = sim_dch.simulate(t_f)  #  78087.336183830  131865.32
            
        SV_dch_df = label_columns(t_dch, SV_dch, an.npoints, sep.npoints, cat.npoints)

        # Obtain tag strings for dataframe columns
        tags = tag_strings(SV_dch_df)
        plot_sim(tags, SV_dch_df, 'Discharging', 0+2*cycle_num, fig, axes)
#        plot_meanPS(SV_dch_df, tags, 'Discharging')
        
        print('Done Discharging\n')
        
        "--------Re-equilibration---------"
        
        if inputs.flag_req == 1:
            
            print('Re-equilibrating...')
            
            # New initial conditions from previous simulation
            SV_0 = SV_dch[-1, :]
            SV_dot_0 = SV_dot_dch[-1, :]
            
            # Set external current
            cat.set_i_ext(0)
            
            # Update problem instance initial conditions
            bat_req = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
            bat_req.external_event_detection = True
            bat_req.algvar = algvar
            
            # Re-initialize simulation object
            sim_req = IDA(bat_req)
            sim_req.atol = atol
            sim_req.rtol = rtol
            sim_req.verbosity = sim_output
            sim_req.make_consistent('IDA_YA_YDP_INIT')
            
            t_req, SV_req, SV_dot_req = sim_req.simulate(t_f)
            
            SV_req_df = label_columns(t_req, SV_req, an.npoints, sep.npoints, cat.npoints)
            
#            plot_sim(tags, SV_req_df, 'Re-Equilibrating', 1, fig, axes)
        
            print('Done re-equilibrating\n')
        else:
            SV_req = SV_dch
            SV_dot_req = SV_dot_dch
            
        "-----------Charging-----------"
        if 'dog' in inputs.ctifile:  
            print('Charging...')
            
            SV_0 = SV_req[-1, :]  
            SV_dot_0 = SV_dot_req[-1, :]  
            
            SV_0[cat.ptr_vec['eps_S8']] = cat.eps_cutoff  
            
            cat.set_i_ext(-cat.i_ext_amp)
            
            # Update problem instance initial conditions
            bat_ch = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
            bat_ch.external_event_detection = True
            bat_ch.algvar = algvar
            
            # Re-initialize simulation object
            sim_ch = IDA(bat_ch)
            sim_ch.atol = atol_ch
            sim_ch.rtol = rtol_ch
            if cycle_num > 0:
                sim_ch.maxh = 0.1
            sim_ch.verbosity = sim_output
            sim_ch.make_consistent('IDA_YA_YDP_INIT')
            
            t_ch, SV_ch, SV_dot_ch = sim_ch.simulate(t_f)
                
            SV_ch_df = label_columns(t_ch, SV_ch, an.npoints, sep.npoints, cat.npoints)
            
            plot_sim(tags, SV_ch_df, 'Charging', 1+inputs.flag_req+2*cycle_num, fig, axes)
            
            plot_meanPS(SV_ch_df, tags, 'Charging')
        SV_ch_df = 0
#        
#        print('Max S_8(e) concentration = ', max(SV_ch[:, 6]))
#        SV_0_cycle = SV_ch[-1, :]
#        SV_dot_0_cycle = SV_ch[-1, :]
        
        print('Done Charging\n')  
    
    exp_data_01C = pd.read_csv(r'0.1C Data.csv', header=None)  
    exp_data_05C = pd.read_csv(r'0.5C Data.csv', header=None)
    exp_data_1C = pd.read_csv(r'1C Data.csv', header=None)
    Bessler = pd.read_csv(r'Bessler Dennis Data.csv', header=None)
    SV_copy = SV_dch_df.copy()
    SV_copy.loc[:, 'Time'] *= -cat.i_ext_amp*inputs.A_cat/3600/(cat.m_S_tot_0)
    "Set up your figure"
    fig = plt.figure(3)
    ax = fig.add_axes([0.2,0.2,0.6,0.75])
    fig.set_size_inches((10.,5.0))
#    ax2 = ax.twinx()
    "Formatting for the figure:"
    fs = 20     #font size for plots
    lw = 3.0    #line width for plots
#    font = plt.matplotlib.font_manager.FontProperties(family='Times New Roman',size=fs-1)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
    color = matplotlib.cm.plasma(inputs.C_rate)
    p1, = plt.plot(SV_copy.loc[:, 'Time'], SV_copy.loc[:, 'Phi_ed1'], 'k-', linewidth=lw)
    p2, = plt.plot(exp_data_01C.iloc[:,0], exp_data_01C.iloc[:,1], 'ro')
    p3, = plt.plot(exp_data_05C.iloc[:,0], exp_data_05C.iloc[:,1], 'co')
    p4, = plt.plot(exp_data_1C.iloc[:,0], exp_data_1C.iloc[:,1], 'ko')
#    p5, = plt.plot(Bessler.iloc[:,0], Bessler.iloc[:,1], 'mo', ms=4)
    plt.xlim((0, 1700))
    plt.xticks([0, 400, 800, 1200, 1600])
    plt.ylim((1.8, 2.6))
    plt.ylabel(r'Cell Voltage $[\mathrm{V}]$', fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
    plt.xlabel(r'Capacity $[\mathrm{Ah} \hspace{0.5} \mathrm{kg}^{-1}_{\mathrm{sulfur}}]$', fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
#        plt.legend(["Discharge", "Charge"])
    th_cat = str(int(inputs.H_cat*1e6))
    file_name_dch = 'dch'+str(inputs.C_rate)+"C_"+th_cat+"um_"+inputs.mech+'.csv'
    SV_dch = SV_dch_df.copy()
    SV_dch.loc[:, 'Time'] *= -cat.i_ext_amp*inputs.A_cat/3600/(cat.m_S_0 + cat.m_S_el)
    SV_dch.to_csv(file_name_dch, index=False, header=True)
    
    t_elapsed = time.time() - t_count
    print('t_cpu=', t_elapsed, '\n')    
    
    return SV_eq_df, SV_dch_df, SV_ch_df, tags 
    
"=============================================================================" 
"===========RESIDUAL CLASSES AND HELPER FUNCTIONS BEYOND THIS POINT==========="
"============================================================================="
from li_s_battery_init import sulfur_obj as sulfur
from li_s_battery_init import Li2S_obj as Li2S
from li_s_battery_init import carbon_obj as carbon
from li_s_battery_init import sulfur_el_s as S_el_s
from li_s_battery_init import Li2S_el_s as L_el_s
from li_s_battery_init import carbon_el_s as C_el_s
from li_s_battery_init import Li2S_tpb
from li_s_battery_init import lithium_obj as lithium
from li_s_battery_init import lithium_el_s as lithium_s
from li_s_battery_init import conductor_obj as conductor
from li_s_battery_init import elyte_obj as elyte
from li_s_battery_functions import set_state, set_geom, set_rxn
from li_s_battery_functions import set_state_sep, set_state_anode
from li_s_battery_functions import dst
from li_s_battery_functions import scale_Diff
from math import pi, exp, tanh

class cc_cycling(Implicit_Problem):    
    def res_fun(t, SV, SV_dot):
        
        res = np.zeros_like(SV)
        ptr = cat.ptr; F = ct.faraday; R = ct.gas_constant; T = inputs.T
        
        """=============================CATHODE============================="""
        """CC BOUNDARY"""
        j = 0; offset = cat.offsets[int(j)]
        i_ext = cat.get_i_ext()
        s2 = set_state(SV, offset, cat.ptr)

        # Set electronic current and ionic current boundary conditions
        i_el_p = i_ext
        i_io_p = 0
        N_io_p = 0  
        
        """=============================CATHODE============================="""
        """INTERIOR NODES"""
        for j in np.arange(1, cat.npoints):
            
            # Set previous outlet fluxes to new inlet fluxes
            i_el_m = i_el_p
            i_io_m = i_io_p
            N_io_m = N_io_p
            s1 = dict(s2)
            
            # Update offset to NEXT node
            offset = cat.offsets[int(j)]
            
            s2 = set_state(SV, offset, cat.ptr)
            
            # Set variables to CURRENT NODE value
            offset = cat.offsets[int(j-1)]
            eps_S8 = max(SV[offset + ptr['eps_S8']], cat.eps_cutoff)
            eps_Li2S = max(SV[offset + ptr['eps_Li2S']], cat.eps_cutoff)
            eps_el = 1 - cat.eps_C_0 - eps_S8 - eps_Li2S  

            # Set states for THIS node            
            carbon.electric_potential = s1['phi_ed']
            elyte.electric_potential = s1['phi_el'] 
            conductor.electric_potential = s1['phi_ed']
            elyte.X = s1['X_k']
#            b = 1e-14  
#            D_scale = b*abs(inputs.C_k_el_0[cat.ptr['iFar']] - s1['C_k'][cat.ptr['iFar']])
            D_scale = scale_Diff(s1['C_k'])
            D_el = (cat.D_el - D_scale)*eps_el**(cat.bruggeman)
#            if i_ext < 0:
#                print(D_scale, cat.D_el, '\n')
            # Current node plus face boundary fluxes
            i_el_p = cat.sigma_eff*(s1['phi_ed'] - s2['phi_ed'])*cat.dyInv
            N_io_p, i_io_p = dst(s1, s2, D_el, cat.dy, cat.dy)
            
            sdot_C = C_el_s.get_net_production_rates(elyte)
            mult = tanh(eps_S8/cat.eps_dropoff)  
            sdot_S8 = S_el_s.get_creation_rates(sulfur) - mult*S_el_s.get_destruction_rates(sulfur)
            sdot_S = S_el_s.get_net_production_rates(elyte)  

            mult = tanh(eps_Li2S/cat.eps_dropoff)  
            sdot_Li2S = L_el_s.get_creation_rates(Li2S) - mult*(L_el_s.get_destruction_rates(Li2S))
            sdot_L = L_el_s.get_net_production_rates(elyte)
            sdot_tpb = Li2S_tpb.get_creation_rates(Li2S) - mult*(Li2S_tpb.get_destruction_rates(Li2S))
            sdot_tpb_el = mult*Li2S_tpb.get_creation_rates(elyte) - Li2S_tpb.get_destruction_rates(elyte)
            
            np_S = inputs.np_S8_init
            np_L = inputs.np_Li2S_init
            
            A_S = 2*pi*np_S*(3*eps_S8/2/np_S/pi)**(2/3)
            A_L = 2*pi*np_L*(3*eps_Li2S/2/np_L/pi)**(2/3)

            r_S = 3*eps_S8/A_S
            r_L = 3*eps_Li2S/A_L
            
            tpb_len = 3*eps_Li2S/(r_L**2)
            A_C = cat.A_C_0 - (pi*np_S*r_S**2) - (pi*np_L*r_L**2)
            
            R_C = sdot_C*A_C
            if eps_S8 < cat.eps_S8_cutoff:
                R_S = 0*sdot_S*A_S
                sw = 0
            else:
                R_S = sdot_S*A_S
                sw = 1
                
#            if eps_Li2S < 1e-5 and i_ext == 0:
#                R_L = 0*sdot_L*A_L + 0*sdot_tpb_el*tpb_len
#                sw2 = 0
#            else:
            R_L = sdot_L*A_L + sdot_tpb_el*tpb_len
            sw2 = 1
    
            i_C = (C_el_s.get_net_production_rates(conductor)*A_C + 
              (Li2S_tpb.get_creation_rates(conductor)*mult - 
               Li2S_tpb.get_destruction_rates(conductor))*tpb_len)
            i_Far = (i_C)*F/cat.dyInv
#            if eps_S8 < 1e-3:
#                print(-i_Far + i_el_m - i_el_p)
            
            # Net rate of formation
            R_net = R_C + R_S + R_L 
            R_net[cat.ptr['iFar']] += (-i_Far + i_el_m - i_el_p)/cat.dy/F
            
            """Calculate change in Sulfur"""                
            res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] 
                                        - sw*sulfur.volume_mole*sdot_S8*A_S)

            """Calculate change in Li2S"""
            res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] 
                                          - sw2*Li2S.volume_mole*(sdot_Li2S*A_L
                                          + sdot_tpb*tpb_len))

            """Calculate change in electrolyte"""
            res[offset + ptr['rho_k_el']] = (SV_dot[offset + ptr['rho_k_el']] - 
            (R_net + (N_io_m - N_io_p)*cat.dyInv)/eps_el 
            + SV[offset + ptr['rho_k_el']]*(- SV_dot[offset + ptr['eps_S8']] 
                                            - SV_dot[offset + ptr['eps_Li2S']])/eps_el)

            """Calculate change in delta-phi double layer"""
            res[offset + ptr['phi_dl']] = (SV_dot[offset + ptr['phi_dl']] - 
            (-i_Far + i_el_m - i_el_p)*cat.dyInv/cat.C_dl/A_C)

            """Algebraic expression for charge neutrality in all phases"""
            res[offset + ptr['phi_ed']] = i_el_m - i_el_p + i_io_m - i_io_p
            
            """Calculate change in S8 nucleation sites"""
            res[offset + ptr['np_S8']] = SV[offset + ptr['np_S8']] - np_S
            
            """Calculate change in Li2S nucleation sites"""
            res[offset + ptr['np_Li2S']] = SV[offset + ptr['np_Li2S']] - np_L

        """=============================CATHODE============================="""
        """SEPARATOR BOUNDARY"""
        
        i_el_m = i_el_p
        i_io_m = i_io_p
        N_io_m = N_io_p
        s1 = dict(s2)
        
        # Shift forward to NEXT node, first separator node (j=0)
        j = 0; offset = sep.offsets[int(j)]
        
        s2 = set_state_sep(SV, offset, sep.ptr)
        
        # Shift back to THIS node, set THIS node outlet conditions
        j = cat.npoints-1; offset = cat.offsets[int(j)]
        
        # Set variables to CURRENT NODE value
#        geom = set_geom(SV, offset, cat.ptr)
        eps_S8 = max(SV[offset + ptr['eps_S8']], cat.eps_cutoff)
        eps_Li2S = max(SV[offset + ptr['eps_Li2S']], cat.eps_cutoff)
        eps_el = 1 - cat.eps_C_0 - eps_S8 - eps_Li2S  
        
        carbon.electric_potential = s1['phi_ed']
        elyte.electric_potential = s1['phi_el']
        conductor.electric_potential = s1['phi_ed']
        elyte.X = s1['X_k']
        
        # Set outlet boundary conditions for THIS node
        i_el_p = 0
        D_scale = scale_Diff(s1['C_k'])  #b*abs(inputs.C_k_el_0[cat.ptr['iFar']] - s1['C_k'][cat.ptr['iFar']])
        D_el = (cat.D_el - D_scale)*eps_el**(cat.bruggeman)
#        if i_ext < 0:
#            print(D_scale, cat.D_el, D_el, t, '\n')

        N_io_p, i_io_p = dst(s1, s2, D_el, cat.dy, sep.dy)

        sdot_C = C_el_s.get_net_production_rates(elyte)
        
        mult = tanh(eps_S8/cat.eps_dropoff)  
        sdot_S8 = S_el_s.get_creation_rates(sulfur) - mult*S_el_s.get_destruction_rates(sulfur)
        sdot_S = S_el_s.get_net_production_rates(elyte)  

        mult = tanh(eps_Li2S/cat.eps_dropoff)  
        sdot_Li2S = L_el_s.get_creation_rates(Li2S) - mult*(L_el_s.get_destruction_rates(Li2S))
        sdot_L = L_el_s.get_net_production_rates(elyte)
        sdot_tpb = Li2S_tpb.get_creation_rates(Li2S) - mult*(Li2S_tpb.get_destruction_rates(Li2S))
        sdot_tpb_el = mult*Li2S_tpb.get_creation_rates(elyte) - Li2S_tpb.get_destruction_rates(elyte)
            
        np_S = inputs.np_S8_init
        np_L = inputs.np_Li2S_init  #SV[offset + ptr['np_Li2S']]
        A_S = 2*pi*np_S*(3*eps_S8/2/np_S/pi)**(2/3)
        A_L = 2*pi*np_L*(3*eps_Li2S/2/np_L/pi)**(2/3)
        
        r_S = 3*eps_S8/A_S
        r_L = 3*eps_Li2S/A_L
                
        tpb_len = 3*eps_Li2S/(r_L**2)
        
        A_C = cat.A_C_0 - (pi*np_S*r_S**2) - (pi*np_L*r_L**2)
        
        R_C = sdot_C*A_C
        if eps_S8 < cat.eps_S8_cutoff:
            R_S = 0*sdot_S*A_S
            sw = 0
        else:
            R_S = sdot_S*A_S
            sw = 1
        
#        if eps_Li2S < 1e-5 and i_ext == 0:
#            R_L = 0*sdot_L*A_L + 0*sdot_tpb_el*tpb_len
#            sw2 = 0
#        else:
        R_L = sdot_L*A_L + sdot_tpb_el*tpb_len
        sw2 = 1

        i_C = (C_el_s.get_net_production_rates(conductor)*A_C + 
              (Li2S_tpb.get_creation_rates(conductor)*mult - 
               Li2S_tpb.get_destruction_rates(conductor))*tpb_len)
        i_Far = (i_C)*F/cat.dyInv
        
        # Net rate of formation
        R_net = R_C + R_S + R_L 
        R_net[cat.ptr['iFar']] += (-i_Far + i_el_m - i_el_p)/cat.dy/F
                                 
        """Calculate change in Sulfur"""                
        res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] 
                                    - sw*sulfur.volume_mole*sdot_S8*A_S)

        """Calculate change in Li2S"""
        res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] 
                                      - sw2*Li2S.volume_mole*(sdot_Li2S*A_L
                                      + sdot_tpb*tpb_len))
                
        """Calculate change in electrolyte"""
        res[offset + ptr['rho_k_el']] = (SV_dot[offset + ptr['rho_k_el']] - 
        (R_net + (N_io_m - N_io_p)*cat.dyInv)/eps_el
        + SV[offset + ptr['rho_k_el']]*(- SV_dot[offset + ptr['eps_S8']] 
                                        - SV_dot[offset + ptr['eps_Li2S']])/eps_el)

        """Calculate change in delta-phi double layer"""
        res[offset + ptr['phi_dl']] = (SV_dot[offset + ptr['phi_dl']] - 
        (-i_Far + i_el_m - i_el_p)*cat.dyInv/cat.C_dl/A_C)

        """Algebraic expression for charge neutrality in all phases"""
        res[offset + ptr['phi_ed']] = i_el_m - i_el_p + i_io_m - i_io_p

        """Calculate change in S8 nucleation sites"""
        res[offset + ptr['np_S8']] = SV[offset + ptr['np_S8']] - np_S
        
        """Calculate change in Li2S nucleation sites"""
        res[offset + ptr['np_Li2S']] = SV[offset + ptr['np_Li2S']] - np_L
        
        """============================SEPARATOR============================"""
        """INTERIOR NODES"""        
        for j in np.arange(1, sep.npoints):
            i_io_m = i_io_p
            N_io_m = N_io_p
            s1 = dict(s2)
            
            # Set NEXT node conditions
            offset = sep.offsets[int(j)]
            s2 = set_state_sep(SV, offset, sep.ptr)
            
            # Shift back to THIS node
            offset = sep.offsets[int(j-1)]
            D_scale = scale_Diff(s1['C_k'])  #b*abs(inputs.C_k_el_0[cat.ptr['iFar']] - s1['C_k'][cat.ptr['iFar']])
            D_el = sep.D_el - D_scale
#            if i_ext < 0:
#                print(D_scale, cat.D_el, D_el, t, '\n')
            
            # THIS node plus face boundary conditions
            N_io_p, i_io_p = dst(s1, s2, D_el, sep.dy, sep.dy)
            
            res[offset + sep.ptr['rho_k_el']] = (SV_dot[offset + sep.ptr['rho_k_el']]
            - (N_io_m - N_io_p)*sep.dyInv/sep.epsilon_el)
            
            res[offset + sep.ptr['phi']] = i_io_m - i_io_p
        
        """============================SEPARATOR============================"""
        """ANODE BOUNDARY"""
        i_io_m = i_io_p
        N_io_m = N_io_p
        s1 = dict(s2)
        
        # Shift forward to NEXT node
        j = 0; offset = an.offsets[int(j)]
        s2 = set_state(SV, offset, an.ptr)
                
        # Shift back to THIS node
        offset = sep.offsets[-1]
        D_scale = scale_Diff(s1['C_k'])  #b*abs(inputs.C_k_el_0[cat.ptr['iFar']] - s1['C_k'][cat.ptr['iFar']])
        D_el = sep.D_el - D_scale
#        if i_ext < 0:
#            print(D_scale, cat.D_el, D_el, t, '\n\n')
        
        # Current node plus face boundary conditions
        N_io_p, i_io_p = dst(s1, s2, D_el, sep.dy, an.dy_el)
#        N_io_p = N_io_p*[1, 1, 1, 0, 0, 0, 0, 0, 0]
#        i_io_p = np.dot(N_io_p, inputs.z_k_el)*F
#        i_io_p = 0
#        N_io_p = 0
        
        res[offset + sep.ptr['rho_k_el']] = (SV_dot[offset + sep.ptr['rho_k_el']]
        - (N_io_m - N_io_p)*sep.dyInv/sep.epsilon_el)
                
        res[offset + sep.ptr['phi']] = i_io_m - i_io_p
        
        """==============================ANODE=============================="""
        """INTERIOR NODES"""
          
        i_io_m = i_io_p
        N_io_m = N_io_p
        i_el_m = 0
        s1 = dict(s2)
        
        j = 0
        offset = an.offsets[int(j)]
        
        i_el_p = i_ext
        i_io_p = 0
        N_io_p = 0
        
        elyte.X = s1['X_k']
        elyte.electric_potential = s1['phi_el']
        lithium.electric_potential = s1['phi_ed']  #SV[offset + an.ptr['phi_ed']]  #
        conductor.electric_potential = s1['phi_ed']  #SV[offset + an.ptr['phi_ed']]  #
        
        sdot_Li = lithium_s.get_net_production_rates(elyte)
        sdot_Far = lithium_s.get_net_production_rates(conductor)
        
        R_net = sdot_Li*an.A_Li*an.dyInv
        i_Far = sdot_Far*an.A_Li*F
        R_net[cat.ptr['iFar']] += (-i_Far + i_el_m - i_el_p)/an.dy/F
        
#        print(sdot_Far, i_ext, '\n')
        
#        res[sep.offsets[-1] + sep.ptr['phi']] = (SV_dot[sep.offsets[-1] + sep.ptr['phi']]
#        + (-i_Far + i_el_m - i_el_p)/an.C_dl/an.A_Li) 
                   
        res[offset + an.ptr['rho_k_el']] = (SV_dot[offset+ an.ptr['rho_k_el']]
        - (R_net + (N_io_m - N_io_p)*an.dyInv)/an.eps_el)  
        
        res[offset + an.ptr['phi_dl']] = (SV_dot[offset + an.ptr['phi_dl']]
        - (-i_Far + i_el_m - i_el_p)/an.C_dl/an.A_Li) 
        
        res[offset + an.ptr['phi_ed']] = SV[offset + an.ptr['phi_ed']] 
        
        """==============================ANODE=============================="""
        """CC BOUNDARY"""

#        print(res, i_ext, t, '\n\n')
#        if i_ext < 0:
#            print(t, '\n\n')
        return res  
    
    "========================================================================="
    
    def state_events(self, t, y, yd, sw):
        
        event1 = np.zeros([cat.npoints])
        event2 = np.zeros([cat.npoints])
        event1 = 1 - y[cat.ptr_vec['eps_S8']]
#        event2 = y[cat.ptr_vec['eps_S8']]
        
        event3 = np.zeros([cat.npoints])
        event4 = np.zeros([cat.npoints])
        event3 = 1 - y[cat.ptr_vec['eps_Li2S']]
#        event4 = y[cat.ptr_vec['eps_Li2S']]
        
        event5 = np.zeros([cat.npoints])
#        event5 = 2.8 - y[cat.ptr_vec['phi_ed']]
        event6 = np.zeros([cat.npoints])
        event6 = y[cat.ptr_vec['phi_ed']] - 1.8
        
        event7 = np.zeros([cat.npoints*elyte.n_species])
        event7 = y[cat.ptr_vec['rho_k_el']]
                
        events = np.concatenate((event1, event2, event3, event4, event5, event6,
                                 event7))

        return events
    
    "========================================================================="
    
    def handle_event(self, solver, event_info):
        state_info = event_info[0]
        
        offset = cat.npoints
        if any(state_info[0:offset]):
            print('Sulfur volume fraction over 1')
            raise TerminateSimulation
        if any(state_info[1*offset:2*offset]):
            print('WARNING: Sulfur volume fraction below 0')
            raise TerminateSimulation
        elif any(state_info[2*offset:3*offset]):
            print('Li2S volume fraction over 1')
            raise TerminateSimulation
        elif any(state_info[3*offset:4*offset]):
            print('Li2S volume fraction below 0')
            raise TerminateSimulation
        elif any(state_info[4*offset:5*offset]):
            print('Cell voltage hit 2.8')
            raise TerminateSimulation
        elif any(state_info[5*offset:6*offset]):
            print('Cell voltage hit 1.8')
            raise TerminateSimulation
        elif any(state_info):
            print('Stop condition')
            raise TerminateSimulation    
    
if __name__ == "__main__":
    SV_eq, SV_dch, SV_ch, tags = main()

"============================================================================="
