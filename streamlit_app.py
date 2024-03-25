from datetime import date
import datetime
import streamlit as st
import numpy as np
import pandas as pd
import altair as alt
import pykep as pk
from matplotlib import pyplot as plt
import numpy as np 
import pandas as pd
import csv
import os
import math
import matplotlib as mpl

# Page title
st.set_page_config(page_title='Mission Analyzer', page_icon='🚀')
st.title('🚀 Mission Analyzer')

with st.expander('About this app'):
  st.markdown('**What can this app do?**')
  st.info("This app computes the Delta-V required to reach a target orbit / body for both first and second burns alongwith the time of flight using Lambert's equations.")
  st.markdown('**How to use the app?**')
  st.warning('To engage with the app, enter initial & final orbit parameters including the mission start and end date from the widget below. As a result, this should generate an updated plot of Delta V calculations based on inputs provided.')



def orbit2orbit_lambert(sc_mass=100, # provide an int or float value in kg
                   sc_thrust=1, # provide an int or float value in N
                   sc_isp=3000, # provide an int or float value in s
                   final_planet_mu = 30): # default value for Ryugu - 30 m3/s2
    
    initial_orbital_elements = [st.session_state.ia, st.session_state.ie, st.session_state.ii*pk.DEG2RAD, st.session_state.iraan*pk.DEG2RAD, st.session_state.iaop*pk.DEG2RAD, st.session_state.ita*pk.DEG2RAD]
    final_orbital_elements = [st.session_state.fa, st.session_state.fe, st.session_state.fi*pk.DEG2RAD, st.session_state.fraan*pk.DEG2RAD, st.session_state.faop*pk.DEG2RAD, st.session_state.fta*pk.DEG2RAD]
    initial_planet = '399'
    final_planet = '399'
    start_date = date.strftime(st.session_state.start, '%Y%m%dT%H%M%S')
    end_date = date.strftime(st.session_state.end, '%Y%m%dT%H%M%S')
    # Initializing spacecraft object
    #sc = pk.sims_flanagan.spacecraft(sc_mass, sc_thrust, sc_isp)
    
    # 1. Getting r, v from keplerian elements for initial orbit in Planet-centric coordinates
    
    # Calculating eccentric anomaly from true anomaly and eccentricity
    true_anomaly_initial = initial_orbital_elements[5] # fetching true anomaly from 
    eccentricity_initial = initial_orbital_elements[1] # fetching the value of orbit eccentricity
    eccentric_anomaly_initial = ((((1 - (eccentricity_initial*eccentricity_initial))**0.5)*(math.sin(true_anomaly_initial)))/(eccentricity_initial + math.cos(true_anomaly_initial)))
    
    # Replacing true anomaly with eccentric anomaly in the list - a, e, i, RAAN, AOP, EA
    initial_orbital_elements[5] = eccentric_anomaly_initial
    
    initial_planet_obj = pk.planet.spice(initial_planet, '0', 'ECLIPJ2000', 'NONE', pk.MU_SUN)
    mu_initial = initial_planet_obj.mu_self

    r_sc, v_sc = pk.par2ic(initial_orbital_elements, mu_initial)

    r_0_sc = np.linalg.norm(r_sc)
    v_0_sc = np.linalg.norm(v_sc)
    print(r_0_sc)
    print(v_0_sc)
    
    # 1. Getting r, v from keplerian elements for initial orbit in Planet-centric coordinates
    
    # Calculating eccentric anomaly from true anomaly and eccentricity
    true_anomaly_final = final_orbital_elements[5]
    eccentricity_final = final_orbital_elements[1]
    eccentric_anomaly_final = ((((1 - (eccentricity_final*eccentricity_final))**0.5)*(math.sin(true_anomaly_final)))/(eccentricity_final + math.cos(true_anomaly_final))) 
    
    # Replacing true anomaly with eccentric anomaly in the list - a, e, i, RAAN, AOP, EA
    final_orbital_elements[5] = eccentric_anomaly_final 
    
    final_planet_obj = pk.planet.spice(final_planet, '0', 'ECLIPJ2000', 'NONE', pk.MU_SUN)
    mu_final = final_planet_mu
    
    r_sc_f, v_sc_f = pk.par2ic(final_orbital_elements, mu_final)
    
    r_1_sc = np.linalg.norm(r_sc_f)
    v_1_sc = np.linalg.norm(v_sc_f)
    print(r_1_sc)
    print(v_1_sc)
    
    #r_sc, v_sc = pk.propagate_lagrangian(r0 = [1,0,0], v0 = [0,1,0], tof = pi/2, mu = 1)


    # START EPOCH 
    start_epoch = pk.epoch_from_iso_string(start_date)  # String format - YYYYMMDDTHHMMSS 
    start_epoch_mjd2000 = start_epoch.mjd2000

     
    # END EPOCH 
    #end_date = "20240110T235858"
    end_epoch = pk.epoch_from_iso_string(end_date)  # String format - YYYYMMDDTHHMMSS 
    end_epoch_mjd2000 = end_epoch.mjd2000

    t1 = start_epoch_mjd2000
    t2 = end_epoch_mjd2000
    dt = (t2 - t1)*86400

    
    r_i, v_i = initial_planet_obj.eph(pk.epoch(start_epoch_mjd2000))
    r_f, v_f = final_planet_obj.eph(pk.epoch(end_epoch_mjd2000))

    r_i = np.array(r_i) + np.array(r_sc)
    v_i = np.array(v_i) + np.array(v_sc)

    r_f = np.array(r_f) + np.array(r_sc_f)
    v_f = np.array(v_f) + np.array(v_sc_f)


    #print(r_l2)
    #print(r_a)

    l = pk.lambert_problem(r1 = r_i, r2 = r_f, tof = dt, mu=pk.MU_SUN, max_revs=0)
    N_max = l.get_Nmax()
    print(N_max)

    print('Computing First Delta-V')
    v10 = l.get_v1()[0]  # Velocity at Departure
    mag_v_l2 = ((v_i[0] ** 2) + (v_i[1] ** 2) + (v_i[2] ** 2)) ** 0.5 # Magnitude of initial Velocity vector from Ephemeris
    mag_v10 = ((l.get_v1()[0][0] ** 2) + (l.get_v1()[0][1] ** 2) + (l.get_v1()[0][2] ** 2)) ** 0.5 # Magnitude of departure Velocity vector
    first_dv = tuple(map(lambda i, j: i - j, v10, v_i))  # First Delta-V vector
    first_dv_magnitude = ((first_dv[0] ** 2) + (first_dv[1] ** 2) + (first_dv[2] ** 2)) ** 0.5 # Magnitude of first Delta-V

    print('Computing Second Delta-V')
    v20 = l.get_v2()[0]
    mag_v_a = ((v_f[0] ** 2) + (v_f[1] ** 2) + (v_f[2] ** 2)) ** 0.5
    mag_v20 = ((l.get_v2()[0][0] ** 2) + (l.get_v2()[0][1] ** 2) + (l.get_v2()[0][2] ** 2)) ** 0.5
    second_dv = tuple(map(lambda i, j: i - j, v20, v_f))
    second_dv_magnitude = ((second_dv[0] ** 2) + (second_dv[1] ** 2) + (second_dv[2] ** 2)) ** 0.5 # Magnitude of second Delta-V

    print('Computing Total dV from initial position to final position')
    total_dv = first_dv_magnitude + second_dv_magnitude

    # Time of flight in days
    tof = t2-t1

    mpl.rcParams['legend.fontsize'] = 10

    # Create the figure and axis
    fig = plt.figure(figsize = (16,8))
    ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    ax1.scatter([0], [0], [0], color=['y'])
    
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax2.scatter([0], [0], [0], color=['y'])
    ax2.view_init(90, 0)
    
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    ax3.scatter([0], [0], [0], color=['y'])
    ax3.view_init(0,0)
    
    for ax in [ax1, ax2, ax3]:
        # Plot the planet orbits
        pk.orbit_plots.plot_planet(initial_planet_obj, t0=t1, color='b', legend=True, units=pk.AU, axes=ax)
        pk.orbit_plots.plot_planet(final_planet_obj, t0=t2, color='#653700', legend=True, units=pk.AU, axes=ax)
        # Plot the Lambert solutions
        axis = pk.orbit_plots.plot_lambert(l, color='g', legend=True, units=pk.AU, axes=ax)
        axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True)
        #axis = pk.orbit_plots.plot_lambert(l, sol=1, color='g', legend=True, units=AU, axes=ax)
        #axis = pk.orbit_plots.plot_lambert(l, sol=2, color='g', legend=True, units=AU, axes=ax)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True)
        
    plt.legend(bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand")
    filename=f"mission_results_{initial_planet}_{final_planet}.png"
    plt.savefig(filename, dpi=500)

    with col2:
      st.subheader('Results')
      with st.expander('Output'):
        st.markdown('Inital Planet')
        st.info(f"Position Vector: {r_i} m")
        st.info(f"Velocity Vector: {v_i} m/sec")
        st.markdown('Final Planet')
        st.info(f"Position Vector: {r_f} m")
        st.info(f"Velocity Vector: {v_f} m/sec")

        st.markdown('First Burn')
        st.info(f"Magnitude: {first_dv_magnitude} m/sec")
        st.info(f"Delta-V vector: {first_dv} m/sec")

        st.markdown('Second Burn')
        st.info(f"Magnitude: {second_dv_magnitude} m/sec")
        st.info(f"Delta-V vector: {second_dv} m/sec")
        st.markdown('Total Delta-V')
        st.success(f"{total_dv} m/sec")
        st.markdown('Time of Flight')
        st.success(f"{tof} days")
      st.pyplot(fig)
      with open(filename, "rb") as file:
        btn = st.download_button(
                label="Download",
                data=file,
                file_name=filename,
                mime="image/png"
              )
    return [initial_planet_obj, final_planet_obj, start_epoch, end_epoch, l]

#sc = [1000, 1, 3000]
col1, col2 = st.columns([0.4,0.6]) 
with col1:
  with st.form(key='my_form'):
    st.subheader('Input Parameters')
    with st.expander('Timeline'):
      st.date_input("Start Date", key="start", value=datetime.date(2026, 1, 1))
      st.date_input("End Date", key="end", value=datetime.date(2027, 1, 1))

    with st.expander('Initial Orbit'):
      st.number_input("A", key="ia")
      st.number_input("E", key="ie")
      st.number_input("I", key="ii",value=23.5)
      st.number_input("RAAN", key="iraan")
      st.number_input("AOP", key="iaop")
      st.number_input("TA", key="ita")

    with st.expander('Final Orbit'):
      st.number_input("A", key="fa")
      st.number_input("E", key="fe")
      st.number_input("I", key="fi", value=23.5)
      st.number_input("RAAN", key="fraan")
      st.number_input("AOP", key="faop")
      st.number_input("TA", key="fta")

    submit_button = st.form_submit_button(label='Submit', type="primary", on_click=orbit2orbit_lambert)

spk_files = "spkfiles/"
print("Loading Spk Files")
pk.util.load_spice_kernel(spk_files + 'de421.bsp')
pk.util.load_spice_kernel(spk_files + 'naif0012.tls')
pk.util.load_spice_kernel(spk_files + '20000399.bsp')
pk.util.load_spice_kernel(spk_files + '20162173.bsp')

# orbit2orbit_lambert(final_orbit, init_orbit, final_planet, init_planet, start_date, end_date)

