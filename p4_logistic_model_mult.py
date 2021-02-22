# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW11 - Problem 4. Logistic model, simulate multiple harvest rates
       K, r and P0 values from previous HW

Created on Sat Feb 13 22:56:55 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt

dt = 1/24  # time step, in days
K = 1.8083  # maximum carrying capacity
r = 0.4477  # low density growth rate
h_rates = [0, r/2, r/4, 3*r/4, 5*r/4]  # harvest rates
labels = ['0', 'r/2', 'r/4', '3r/4', '5r/4']  # corresponding labels for plotting
P0 = 0.0498  # initial population
steps = 50

t = np.linspace(0, steps, int(steps/dt)+1)  # time vector
print(t)

# Make predictions using the analytical model and numpy array
Pexp = (K * P0 * np.exp(r*t))/(K + P0 * (np.exp(r*t) - 1))

# A list to contain all the simulations with a different value of 'h'
Pharvest = []

for h in h_rates:
    # initialize the population vectors
    Pharv = np.zeros(len(t))
    Pharv[0] = P0
    
    for i in range(1, len(t)):
        Pharv[i] = Pharv[i-1] + r * Pharv[i-1] * (1 - Pharv[i-1]/K)*dt - (h * Pharv[i-1] * dt)
    Pharvest.append(Pharv)

plt.figure(0)
plt.title('HW11 Algae model with harvest')
plt.xlabel('Time (days)')
plt.ylabel('Population')
# plt.plot(t, Pexp, 'r.-', label='Analytic (eq.)')  # Analytical solution
# Plot every simulation
for Ph, lbl in zip(Pharvest, labels):
    plt.plot(t, Ph, label='h=' + lbl)

plt.legend(loc='best')
plt.savefig('p4_logistic_model_%ddays_mult.png' % steps, dpi=300, bbox_inches='tight')
plt.show()
