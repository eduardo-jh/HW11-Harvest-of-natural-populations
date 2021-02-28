# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW11 - Problem 4. Logistic model, harvest when concentration reaches K/2
       K, r and P0 values from previous HW

Created on Sat Feb 13 22:56:55 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt

dt = 1/24  # time step, in days
K = 1.8083  # maximum carrying capacity
r = 0.4477  # low density growth rate
hr = r/2  # harvest rate
P0 = 0.0498  # initial population
steps = 50

t = np.linspace(0, steps, int(steps/dt)+1)  # time vector
print(t)

# initialize the population vectors
P = np.zeros(len(t))
Pharv = np.zeros(len(t))
P[0], Pharv[0] = P0, P0

for i in range(1, len(t)):
    P[i] = P[i-1] + r * P[i-1] * (1 - P[i-1]/K)*dt
    h = (hr if Pharv[i-1] > K/2 else 0)  # harvest when concentration reaches K/2
    Pharv[i] = Pharv[i-1] + r * Pharv[i-1] * (1 - Pharv[i-1]/K)*dt - (h * Pharv[i-1] * dt)

# Make predictions using the analytical model and numpy array
Pexp = (K * P0 * np.exp(r*t))/(K + P0 * (np.exp(r*t) - 1))

plt.plot(t, P, 'b-', t, Pharv, 'g-', t, Pexp, 'r-')
plt.legend(['Euler (dt=%.1f h)' % (dt*24), 'Euler (harvest h=%.4f)' % hr, 'Analytic (eq.)'], loc='best')
plt.title('HW11 Algae model with harvest')
plt.xlabel('Time (days)')
plt.ylabel('Population')
# plt.savefig('p4_logistic_model_%ddays_0.5K.png' % steps, dpi=300, bbox_inches='tight')
plt.show()
