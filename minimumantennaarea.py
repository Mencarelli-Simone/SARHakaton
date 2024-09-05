import numpy as np
from numpy import arctan, tan, sin, cos, vectorize, arctan2, arccos, pi
from design_functions import *
from spherical_earth_geometry_radar import *

# %% User input
freq = 10e9  # 5e9
La = 6  # antenna length

# incidence angle
eta = 35 * np.pi / 180

# altitude
h = 500e3

# dutycycle
dtc = 20 / 100

# ground range resolution
rrg = 3

# losses, noise figure, efficiency (i.e.  power budget)
Loss = 10  # dB

# NESZ level goal
NESZ = -20  # dB

print('power budget assumptions:')
print('Loss + Nfigure + efficiency: {} dB'.format(Loss))
print('NESZ: {} dB'.format(NESZ))
print('operating frequency: {:.2f} GHz'.format(freq / 1e9))
print('incidence angle: {:.2f} deg'.format(eta * 180 / np.pi))
print('Antenna length: {:.2f} m'.format(La))

# %%
# nominal dopplere
bd = nominal_doppler_bandwidth(La, eta, 3e8 / freq, orbital_speed(h), h=500e3)
# doppler oversampling
osd = 1.1
PRF = bd * osd
print('Nominal Doppler bandwidth: {:.2f} Hz'.format(bd))
print('PRF: {:.2f} Hz'.format(PRF))
# swath
rs, _ = range_from_theta(eta * 180 / np.pi, h=h)
rank = int(2 * rs * PRF / 3e8)
print('rank: {}'.format(rank))
rne = rank * 3e8 / (2 * PRF) + dtc / 2 * 3e8 / (2 * PRF)
rfe = (rank + 1) * 3e8 / (2 * PRF) - dtc / 2 * 3e8 / (2 * PRF)
rg, et = range_slant_to_ground(np.array([rne, rfe]), h=h)
wg = rg[1] - rg[0]
print('Swath: {:.2f} km'.format(wg / 1e3))
# beamwidth
loka = incidence_angle_to_looking_angle(et, h=h)
dang = loka[1] - loka[0]
Wa = 3e8 / (freq * dang)
print('Beamwidth: {:.2f} deg'.format(dang * 180 / np.pi))
print('Antenna width: {:.2f} m'.format(Wa))

# range bandwidth
Bn = find_bandwidth(La, eta * 180 / np.pi, La / 2 * rrg)
print('Range bandwidth: {:.2f} MHz'.format(Bn / 1e6))

# NESZ level
# gain of the antenna
G = 4 * pi * (La * Wa) / (3e8 / freq) ** 2
print('Antenna gain: {:.2f} dB'.format(10 * np.log10(G)))
vg = ground_speed(np.mean(rg), orbital_speed(h), h=h)
# boltzmann constant
k_boltz = 1.380649E-23  # J/K
# antenna temperature
Tant = 300
Pavg = rs ** 3 * 256 * pi ** 3 * Bn * sin(eta) * vg * k_boltz * Tant * 10 ** (Loss / 10) / (10 ** (NESZ / 10) *
                                                                                            G ** 2 * (
                                                                                                        3e8 / freq) ** 3 * 3e8)
print('Required Average power: {:.2f} W'.format(Pavg))
