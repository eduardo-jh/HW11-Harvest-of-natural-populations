#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW11 - Problem 5. Multiple fission algae model, with harvest simulation

Created on Sun Feb 14 00:40:51 2021
@author: eduardo
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

number_of_days = 30  # days
kN = 3e-5  # Nitrogen, mole/L
kP = 1e-5  # Phophorous, mole/L
kCO2 = 2.5e-5  # Carbon dioxide, mole/L
time_step = 1 # time step, hours
C0 = 0.05  # Initial concentration, mg/L
max_light = 800  # W/m2
k = 0.002  # Light mu multiplier, W/m2
mu_night = -0.1  # 1/day
max_conc = 2 #  mg/L
N0 = 0.1  # Nitrogen initial, moles/L
P0 = 0.01  # Phophorous initial, moles/L
CO20 = 1  # Carbon dioxide initial, moles/L

r = 0.4477  # low density growth rate
h = r/2  # harvest rate

hpd = 24  # hours per day
N_moles = 16
P_moles = 1
CO2_moles = 124

dt = time_step * (1/hpd)  # convert time step into days
time_steps = int(number_of_days / dt)  # get time steps in the number of days
# This is grams of algae per mole of P, or 16 moles of N
algae_conc = 12*106 + 1*263 + 16*110 + 14*16 + 31  # total g per mole

time = np.linspace(0, (time_steps*dt), time_steps+1)  # time vector

C = np.ones(time_steps+1) * C0
N = np.ones(time_steps+1) * N0
P = np.ones(time_steps+1) * P0
CO2 = np.ones(time_steps+1) * CO20
mu = np.ones(time_steps+1) * mu_night

algae_yield = np.zeros(time_steps+1)
I0 = np.zeros(time_steps+1)

for i in range(1, time_steps+1):
    I0[i] = max_light * np.sin(2 * np.pi * (i-6)/hpd)
    I0_effective = I0[i] * (max_conc - C[i-1]) / max_conc
    
    mu[i] = I0_effective * k * N[i-1]/(kN + N[i-1]) * P[i-1]/(kP + P[i-1]) * CO2[i-1]/(kCO2 + CO2[i-1])
    
    mu[i] = (mu_night if mu[i] < 0 else mu[i])  # use night coefficient when negative values come out

    # Harvest the algae after reaching half of maximum concentration 'K'
    hr = (h if C[i-1] > (max_conc/2) else 0) # choose a different harvest rate before and after K/2
    deltaC = C[i-1] * (mu[i] - hr) * dt
    C[i] = C[i-1] + deltaC
    algae_yield[i] = C[i-1] * h * dt #harvested algae
    
    if mu[i] > 0:
        N[i] = N[i-1] - deltaC / algae_conc * N_moles
        P[i] = P[i-1] - deltaC / algae_conc * P_moles
        CO2[i] = CO2[i-1] - deltaC / algae_conc * CO2_moles
    else:
        N[i] = N[i-1]
        P[i] = P[i-1]
        CO2[i] = CO2[i-1]

growth_script = pd.DataFrame(columns = ['Time', 'I0', 'mu', 'Conc', 'N', 'P', 'CO2'])
growth_script['Time'] = time
growth_script['I0'] = I0
growth_script['mu'] = mu
growth_script['Conc'] = C
growth_script['N'] = N
growth_script['P'] = P
growth_script['CO2'] = CO2

plt.figure(0)
plt.plot(growth_script['Time'], growth_script['Conc'], 'k-', label='Conc')
plt.plot(growth_script['Time'], growth_script['N'], 'b--', label='N')
plt.plot(growth_script['Time'], growth_script['P'], 'r-.', label='P')
plt.plot(growth_script['Time'], growth_script['CO2'], 'm:', label='$CO_2$')
plt.plot(growth_script['Time'], growth_script['mu'], 'c-', label='$\mu$')
plt.legend(loc='best', fancybox=True)
plt.xlabel('Time [days]')
plt.ylabel('Concentration $K,P,CO_2$ [mole/L], C [mg/L]')
# plt.savefig('p5_algae_harvest_%ddays.png' % number_of_days, dpi=300, bbox_inches='tight')
plt.show()
