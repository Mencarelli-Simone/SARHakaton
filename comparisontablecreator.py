import numpy as np
from numpy import pi, sin
from design_functions import *
from spherical_earth_geometry_radar import *


# Design function
def radar_design(freq, La, eta, h, dtc, rrg, Loss, NESZ):
    """
    Perform radar design calculations based on the provided requirements.

    Parameters:
        freq (float): Carrier frequency (Hz).
        La (float): Antenna length (m).
        eta (float): Incidence angle (rad).
        h (float): Altitude (m).
        dtc (float): Duty cycle (fraction, e.g., 0.2 for 20%).
        rrg (float): Ground range resolution (m).
        Loss (float): Loss + noise figure + efficiency (dB).
        NESZ (float): Noise Equivalent Sigma Zero (dB).

    Returns:
        dict: A dictionary containing the calculated design outputs.
    """
    # Perform calculations
    bd = nominal_doppler_bandwidth(La, eta, 3e8 / freq, orbital_speed(h), h=h)
    it = integration_time(La, eta, 3e8 / freq, orbital_speed(h), h=h)  # Integration time
    PRF = bd * 1.1  # Doppler oversampling factor
    rs, _ = range_from_theta(eta * 180 / pi, h=h)
    rank = int(2 * rs * PRF / 3e8)
    rne = rank * 3e8 / (2 * PRF) + dtc / 2 * 3e8 / (2 * PRF)
    rfe = (rank + 1) * 3e8 / (2 * PRF) - dtc / 2 * 3e8 / (2 * PRF)
    rg, et = range_slant_to_ground(np.array([rne, rfe]), h=h)
    wg = rg[1] - rg[0]
    loka = incidence_angle_to_looking_angle(et, h=h)
    dang = loka[1] - loka[0]
    Wa = 3e8 / (freq * dang)
    Bn = find_bandwidth(La, eta * 180 / pi, La / 2 * rrg)
    G = 4 * pi * (La * Wa) / (3e8 / freq) ** 2
    vg = ground_speed(np.mean(rg), orbital_speed(h), h=h)
    k_boltz = 1.380649E-23  # Boltzmann constant (J/K)
    Tant = 300  # Antenna temperature (K)
    Pavg = rs ** 3 * 256 * pi ** 3 * Bn * sin(eta) * vg * k_boltz * Tant * 10 ** (Loss / 10) / (
            10 ** (NESZ / 10) * G ** 2 * (3e8 / freq) ** 3 * 3e8
    )

    # Return design outputs
    return {
        "$W_g$ (km)": wg / 1e3,
        "$B_D$ (Hz)": bd,
        "PRF (Hz)": PRF,
        "$T_i$ (s)": it,  # Integration time
        "$B_n$ (MHz)": Bn / 1e6,
        "$\\Theta_{el}$ (deg)": dang * 180 / pi,
        "$W_A$ (m)": Wa,
        "G (dB)": 10 * np.log10(G),
        "$P_{av}$ (W)": Pavg,
    }


# Function to check for unchanged parameters and collapse them
def collapse_parameters(data, toggle):
    """
    Collapse unchanged parameters across designs.

    Parameters:
        data (dict): Dictionary containing parameters and their values across designs.
        toggle (bool): Whether to collapse unchanged parameters.

    Returns:
        dict: Modified dictionary with collapsed parameters if enabled.
    """
    if not toggle:
        return data  # Return unchanged if toggle is off

    for key, values in data.items():
        if all(v == values[0] for v in values):  # Check if all values are identical
            data[key] = [values[0]] + ["-"] * (len(values) - 1)
    return data


# Function to write results to a LaTeX table
def write_comparison_table(requirements, designs, filename="comparison_results.tex"):
    with open(filename, "w") as f:
        # Start table
        f.write("\\begin{table}[h!]\n")
        f.write("\\centering\n")
        f.write("\\caption{Comparison of Radar Designs}\n")
        f.write("\\begin{tabular}{|l|" + "c|" * len(designs) + "}\n")
        f.write("\\hline\n")

        # Header with design indices
        header = ["\\textbf{Parameter}"] + [f"\\textbf{{Design {i + 1}}}" for i in range(len(designs))]
        f.write(" & ".join(header) + " \\\\\n")
        f.write("\\hline\n")

        # Requirements section
        f.write("\\multicolumn{" + str(len(designs) + 1) + "}{|c|}{\\textbf{Requirements}} \\\\\n")
        f.write("\\hline\n")
        for param, values in requirements.items():
            row = [param] + values
            f.write(" & ".join(map(str, row)) + " \\\\\n")

        # Design outputs section
        f.write("\\hline\n")
        f.write("\\multicolumn{" + str(len(designs) + 1) + "}{|c|}{\\textbf{Design Outputs}} \\\\\n")
        f.write("\\hline\n")
        for param in designs[0].keys():
            row = [param] + [f"{design[param]:.2f}" for design in designs]
            f.write(" & ".join(row) + " \\\\\n")

        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")


# Main section
if __name__ == "__main__":
    # Toggle for collapsing unchanged parameters
    collapse_unchanged = True

    # Input variations
    frequencies = [10e9, 10e9, 5.4e9, 35.8e9]  # Carrier frequencies (Hz)
    antenna_lengths = [2, 4.8, 15, 5]  # Antenna lengths (m)
    assert len(frequencies) == len(antenna_lengths), "Frequencies and antenna lengths must have the same length."

    eta = 30 * pi / 180  # Incidence angle (rad)
    h = 500e3  # Altitude (m)
    dtc = 20 / 100  # Duty cycle
    rrg = 2  # Ground range resolution (m)
    Loss = 10  # Loss + Noise Figure + Efficiency (dB)
    NESZ = -27  # NESZ (dB)

    # Store requirements for each design
    requirements = {
        "$f_c$ (GHz)": [freq / 1e9 for freq in frequencies],
        "$L_A$ (m)": antenna_lengths,
        "$\\eta_0$ (deg)": [eta * 180 / pi] * len(frequencies),
        "h (km)": [h / 1e3] * len(frequencies),
        "Duty Cycle (\\%)": [dtc * 100] * len(frequencies),
        "$\\delta_{r_g}$ (m)": [rrg] * len(frequencies),
        "NESZ (dB)": [NESZ] * len(frequencies),
        "L + N (dB)": [Loss] * len(frequencies),
    }

    # Collapse unchanged parameters if the toggle is on
    requirements = collapse_parameters(requirements, collapse_unchanged)

    # Perform designs
    designs = []
    for i in range(len(frequencies)):
        design = radar_design(frequencies[i], antenna_lengths[i], eta, h, dtc, rrg, Loss, NESZ)
        designs.append(design)

    # Write LaTeX table
    write_comparison_table(requirements, designs)
    print("Comparison table written to 'comparison_results.tex'")
