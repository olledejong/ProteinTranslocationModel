from math import log

# -------------------- Simulation parameters

num_cycles = 6  # number of cycles to include in the multiple-cycles plot
num_datapoints = 100  # the desired number of datapoints that is solved for within the time-axis
nuc_div_tp = 91  # simulated point at which nuc division takes place (based on volume drop-off time-point of nuc vol)

# -------------------- Initial conditions

cp0 = 200  # initial cytosolic protein abundance
np0 = 10  # initial nuclear protein abundance

# -------------------- Rates

ks = 0.25  # synthesis rate of protein in cytosol
kIn = log(2)  # rate of translocation into nucleus
kOut = log(2) / 10  # rate of translocation out of nucleus
kd = log(2) / 35  # degradation rate for protein
kIn_mp = 1.5  # the higher this is, the higher the maximum of the nuclear import rate curve
kOut_mp = 4  # the higher this is, the higher the maximum of the nuclear export rate curve
