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
  
st.subheader('Input Parameters')



def orbit2orbit_lambert(initial_orbital_elements, # [a, e, i, RAAN, AOP, TA]
                   final_orbital_elements, # [a, e, i, RAAN, AOP, TA]
                   initial_planet, # provide a string of the planet's SPK ID. For example for Earth: '399'
                   final_planet, # provide a string of the planet's SPK ID. For example for Moon: '301'
                   start_date, # provide a string of date-time in the format - 'YYYYMMDDTHHMMSS'
                   end_date, # provide a string of date-time in the format - 'YYYYMMDDTHHMMSS'
                   sc_mass=100, # provide an int or float value in kg
                   sc_thrust=1, # provide an int or float value in N
                   sc_isp=3000, # provide an int or float value in s
                   final_planet_mu = 30): # default value for Ryugu - 30 m3/s2
    
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


    """Dates for the simulation"""

    # START EPOCH 
    start_epoch = pk.epoch_from_iso_string(start_date)  # String format - YYYYMMDDTHHMMSS 
    start_epoch_mjd2000 = start_epoch.mjd2000

     
    # END EPOCH 
    #end_date = "20240110T235858"
    end_epoch = pk.epoch_from_iso_string(end_date)  # String format - YYYYMMDDTHHMMSS 
    end_epoch_mjd2000 = end_epoch.mjd2000


    """ Boundary Conditions for Lambert's problem """
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


    """ Delta-V calculations """

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

    print("Computation done!")
    print('-------------------------------------------------------------')
    print("INPUTS GIVEN")
    print('-------------------------------------------------------------')
    print('Start date:', start_epoch)
    print('End date:', end_epoch)
    print('-------------------------------------------------------------')
    print("OUTPUTS")
    print('-------------------------------------------------------------')
    print('position vector at initial planet:', r_i, " m ")
    print('velocity vector at initial planet:', v_i, " m/s ")
    print('position vector at final planet:', r_f, " m ")
    print('velocity vector at final planet', v_f, " m/s ")
    print('First burn delta-V vector:', first_dv, " m/s ")
    print('Magnitude of First burn delta-V:', first_dv_magnitude, " m/s ")
    print('Second burn delta-V vector:', second_dv, " m/s ")
    print('Magnitude of Second burn delta-V:', second_dv_magnitude, " m/s ")
    print('----------------------FINAL VALUES---------------------------')
    print('Total delta-V magnitude:', total_dv, " m/s ")
    print('Time of Flight:', tof, " days ")
    print('-------------------------------------------------------------')

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
    plt.savefig('mission_results'+ '_' + initial_planet + '_' + final_planet, dpi=500)
    st.pyplot(fig)
    return [initial_planet_obj, final_planet_obj, start_epoch, end_epoch, l]

init_orbit = [7000000, 0, 23.5*pk.DEG2RAD, 0*pk.DEG2RAD, 0*pk.DEG2RAD, 0*pk.DEG2RAD]
final_orbit = [10000, 0, 23.5*pk.DEG2RAD, 0*pk.DEG2RAD, 0*pk.DEG2RAD, 0*pk.DEG2RAD]
init_planet = '399'
final_planet = '399'
start_date = "20261001T000000"
end_date = "20270528T235852"
#sc = [1000, 1, 3000]

spk_files = "spkfiles/"
print("Loading Spk Files")
pk.util.load_spice_kernel(spk_files + 'de421.bsp')
pk.util.load_spice_kernel(spk_files + 'naif0012.tls')
pk.util.load_spice_kernel(spk_files + '20000399.bsp')
pk.util.load_spice_kernel(spk_files + '20162173.bsp')

orbit2orbit_lambert(final_orbit, init_orbit, final_planet, init_planet, start_date, end_date)

