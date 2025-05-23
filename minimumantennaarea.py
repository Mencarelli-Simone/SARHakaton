import numpy as np
from numpy import arctan, tan, sin, cos, vectorize, arctan2, arccos, pi
from design_functions import *
from spherical_earth_geometry_radar import *

# %% User input
freq = 3e9  # 5e9
La = 1 # antenna length

# incidence angle
eta = 40 * np.pi / 180

# altitude
h = 20e3

# speed
# replace vs everywhere
# vs = orbital_speed(h)
vs = 800 * 1000/3600

# dutycycle
dtc = 20 / 100

# ground range resolution
rrg = 0.5

# losses, noise figure, efficiency (i.e.  power budget)
Loss = 10  # dB

# NESZ level goal
NESZ = -27# dB

print('power budget assumptions:')
print('Loss + Nfigure + efficiency: {} dB'.format(Loss))
print('NESZ: {} dB'.format(NESZ))
print('operating frequency: {:.2f} GHz'.format(freq / 1e9))
print('incidence angle: {:.2f} deg'.format(eta * 180 / np.pi))
print('Antenna length: {:.2f} m'.format(La))


# %%
# nominal dopplere
bd = nominal_doppler_bandwidth(La, eta, 3e8 / freq, vs, h=500e3)
it = integration_time(La, eta, 3e8 / freq, vs, h=500e3)
# doppler oversampling
osd = 1.1
PRF = bd * osd
print('Nominal Doppler bandwidth: {:.2f} Hz'.format(bd))
print('Integration time: {:.2f} s'.format(it))
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
vg = ground_speed(np.mean(rg), vs, h=h)
# boltzmann constant
k_boltz = 1.380649E-23  # J/K
# antenna temperature
Tant = 300
Pavg = rs ** 3 * 256 * pi ** 3 * Bn * sin(eta) * vg * k_boltz * Tant * 10 ** (Loss / 10) / (10 ** (NESZ / 10) *
                                                                                            G ** 2 * (
                                                                                                        3e8 / freq) ** 3 * 3e8)
print('Required Average power: {:.2f} W'.format(Pavg))

#%% latex table generation
# Collect requirements and design outputs
requirements = {
    '$f_c$ (GHz)': freq / 1e9,
    '$\\eta_0$ (deg)': eta * 180 / np.pi,
    '$L_A$ (m)': La,
    'h (km)': h / 1e3,
    'Duty Cycle (\\%)': dtc * 100,
    '$\\delta_{r_g}$ (m)': rrg,
    'NESZ (dB)': NESZ,
    'L + N (dB)': Loss
}

design_outputs = {
    '$W_g$ (km)': wg / 1e3,
    '$B_D$ (Hz)': bd,
    'PRF (Hz)': PRF,
    '$B_n$ (MHz)': Bn / 1e6,
    '$\Theta_{el}$ (deg)': dang * 180 / np.pi,
    '$W_A$ (m)': Wa,
    'G (dB)': 10 * np.log10(G),
    '$P_{av}$ (W)': Pavg
}


# Function to write results to a LaTeX table
def write_latex_table(reqs, outputs, filename="design_results.tex"):
    with open(filename, "w") as f:
        f.write("\\begin{table}[h!]\n")
        f.write("\\centering\n")
        f.write("\\caption{Design Requirements and Outputs}\n")
        f.write("\\begin{tabular}{|l|c|}\n")
        f.write("\\hline\n")
        f.write("\\textbf{Requirements} & \\textbf{Value} \\\\\n")
        f.write("\\hline\n")

        # Write requirements
        for key, value in reqs.items():
            if isinstance(value, float):
                f.write(f"{key} & {value:.2f} \\\\\n")
            else:
                f.write(f"{key} & {value} \\\\\n")

        f.write("\\hline\n")
        f.write("\\textbf{Design Outputs} & \\textbf{Value} \\\\\n")
        f.write("\\hline\n")

        # Write design outputs
        for key, value in outputs.items():
            if isinstance(value, float):
                f.write(f"{key} & {value:.2f} \\\\\n")
            else:
                f.write(f"{key} & {value} \\\\\n")

        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")


# Call the function to generate the LaTeX table
write_latex_table(requirements, design_outputs)

print("LaTeX table written to 'design_results.tex'")
