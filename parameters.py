from math import log

# -------------------- Simulation parameters

num_cycles = 6  # number of cycles to run
num_datapoints = 100  # the desired number of datapoints that is solved for within the time-axis
nuc_div_tp = 91  # simulated point at which nuc division takes place (based on volume drop-off time-point of nuc vol)
normalize = False  # whether to normalize the compared traces or not

# -------------------- Initial conditions

# Are not of great importance when running multiple consecutive simulations, since the model (in this form) will
# migrate to a stable oscillatory pattern. Once you perform only a single simulation, then they do in fact influence the
# model output.

cp0 = 200  # initial cytosolic protein abundance
np0 = 10  # initial nuclear protein abundance

# -------------------- Rates

ks = 5  # synthesis rate
kIn = log(2)  # rate of translocation into nucleus
kOut = log(2) / 10  # rate of translocation out of nucleus
kd = log(2) / 35  # degradation rate for protein
kIn_mp = 1.5  # the higher this is, the higher the maximum of the nuclear import rate curve
kOut_mp = 4  # the higher this is, the higher the maximum of the nuclear export rate curve
