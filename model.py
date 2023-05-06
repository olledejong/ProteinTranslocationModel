import sys
import plotting
import numpy as np
import pandas as pd
from math import log
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline, interp1d, UnivariateSpline

peak_of_nuc_vol = 81  # the latest local high of the nuclear volume
nuc_div_tp = 91  # simulated point at which nuclear division takes place


########################
### Global variables ###
########################
cell_vols, nuc_vols, nuc_surface_areas, t_range = [], [], [], None
cell_vols: list
nuc_vols: list
nuc_surface_areas: list
a_func, cv_func, nv_func = None, None, None


#########################
### Data preparations ###
#########################
def load_and_adjust_data():
    """
    Loads data from the averages file, calls the alter_data function and returns the outcome.
    :return:
    """
    file = "./averages.xlsx"
    averages = pd.read_excel(file)

    return alter_data(averages)


def alter_data(averages):
    """
    Takes the by the script exported averages and alters these based on various assumptions:
      - Nuclear division takes place approximately 10 time-steps before the end (which represents whole-cell splitting)
      - Nuclear division duration is very short, and thus assumed to be instant (split off at one time-point)
      - Nuclear volume and surf area stays constant from local high up until division point
    :param averages:
    :return:
    """
    global cell_vols, nuc_vols, nuc_surface_areas

    # nuclear volume over time (altered to simulate instant nuclear division)
    nv = averages.nuc_volumes.values
    nv[peak_of_nuc_vol:nuc_div_tp] = nv[peak_of_nuc_vol]
    nv[nuc_div_tp:] = [nv[-1]] * len(nv[nuc_div_tp:])

    # nuclear surface area over time (altered such that it simulates a proportional decrease compared to nuclear volume)
    nsa = averages.nuc_surface_areas.values
    nsa[peak_of_nuc_vol:nuc_div_tp] = nsa[peak_of_nuc_vol]
    nsa[nuc_div_tp:] = [nsa[-1]] * len(nsa[nuc_div_tp:])

    # whole-cell volumes
    cv = averages.cell_volumes.values

    cell_vols = cv
    nuc_vols = nv
    nuc_surface_areas = nsa


def define_area_vol_functions():
    global a_func, cv_func, nv_func, t_range
    # interpolated vol/area functions for time dependent retrieval of values
    t_range = np.arange(len(nuc_surface_areas))
    a_func = interp1d(t_range, nuc_surface_areas, kind='linear', bounds_error=False)
    cv_func = interp1d(t_range, cell_vols, kind='linear', bounds_error=False)
    nv_func = interp1d(t_range, nuc_vols, kind='linear', bounds_error=False)



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


#################
### The model ###
#################
def dp_dt(y, t, k_d, k_s, k_in, k_out):
    C, N = y  # cellular and nuclear abundance of protein
    A = a_func(t)  # nuclear surface area at t
    Vc = cv_func(t)  # whole cell volume at t
    Vn = nv_func(t)  # nuclear volume at t

    # equations that describe the change in cellular and nuclear abundances
    # dC/dt: (1) synthesis of protein (2) net transfer into the cyt (prop to nuc surf area) (3) deg. of prot in cyt
    # dN/dt: (1) net transfer into the nucleus (proportional to nuc surf. area) (2) degradation of prot in nucleus
    dC_dt = k_s * 20 + A * (-k_in * C / (Vc - Vn) + k_out * N / Vn) - k_d * C
    dN_dt = A * (k_in * C / (Vc - Vn) - k_out * N / Vn) - k_d * N

    return [dC_dt, dN_dt]


###########################
### Running simulations ###
###########################
def main():
    load_and_adjust_data()
    define_area_vol_functions()

    num_cycles = 5  # number of cycles to include in the multiple-cycles plot
    num_datapoints = 200  # the desired number of datapoints that is solved for within the time-axis
    vol_loss_frac = 0.28125901660197855  # fraction of total volume that on average is lost to the daughter bud
    time_scalar = 0.7135  # scalar based on average duration of cycle to scale back to real minute axis

    ### Model parameters ###

    kc = 0.25  # synthesis rate of protein in cytoplasm
    kd = log(2) / 35  # degradation rate for protein
    kIn = log(2) / 10 / np.average(nuc_surface_areas)  # rate of translocation into nucleus
    kOut = log(2) / 10 / np.average(nuc_surface_areas)  # rate of translocation out of nucleus
    # log(2)/10 is the rate of translocation, this is scaled by dividing it by the average nuclear surface

    # initial conditions
    cp0 = 1000
    np0 = 60

    ### Running simulations ###
    mult_cycles_cyt, mult_cycles_nuc = [], []
    for i in range(num_cycles):
        # generate split t_span to be able to simulate nuclear division event at 'nuc_div_tp'
        tspan_whole = np.linspace(0, 99, num_datapoints)
        closest = find_nearest(tspan_whole, nuc_div_tp)
        tspan_before, tspan_after = np.split(tspan_whole, np.where(tspan_whole == closest)[0] + 1)

        # solve the ode up until the event
        sols_before = odeint(dp_dt, [cp0, np0], tspan_before, args=(kd, kc, kIn, kOut))
        sols_bf_ev_cyt = sols_before[:, 0]
        sols_bf_ev_nuc = sols_before[:, 1]

        # simulate event (reduction of nuclear abundance proportional to loss in nuc vol)
        perc_to_remove = (np.amax(nuc_vols) - nuc_vols[-1]) / np.amax(nuc_vols)
        sols_after_event = np.array([
            sols_bf_ev_cyt[-1],  # nothing happens to the cellular abundance at nuclear split
            sols_bf_ev_nuc[-1] - (sols_bf_ev_nuc[-1] * perc_to_remove)  # remove the calc. percentage from the ab.
        ])

        # simulate part after the event
        sols_after = odeint(dp_dt, sols_after_event, tspan_after, args=(kd, kc, kIn, kOut))
        sols_after_cyt = sols_after[:, 0][~np.isnan(sols_after[:, 0])]
        sols_after_nuc = sols_after[:, 1][~np.isnan(sols_after[:, 1])]

        # simulate the bud separation event by reducing the cytoplasmic abundance in proportion to the volume loss
        # percentage (0.281..) determined using the volume analysis pipeline
        sols_after_cyt[-1] = (1 - vol_loss_frac) * sols_after_cyt[-1]

        # nans are produced because the outcome of the n_func / cv_func / nv_func is unknown at the last few timepoints
        num_nans = np.count_nonzero(np.isnan(sols_after[:, 0]))

        final_tspan = np.concatenate((tspan_before, tspan_after))  # merge the two split t-spans
        # according to the volume analysis script, the average duration of one cycle is 71.35 minutes.
        # let's modify the time-axis using this knowledge (for purpose of plotting on a real time axis)
        final_tspan = final_tspan * time_scalar

        final_cyt_ab = np.concatenate((sols_before[:, 0], sols_after_cyt))
        final_nuc_ab = np.concatenate((sols_before[:, 1], sols_after_nuc))
        mult_cycles_cyt.extend(final_cyt_ab)
        mult_cycles_nuc.extend(final_nuc_ab)

        # set the initial condition for the next loop to the final values
        cp0, np0 = final_cyt_ab[-1], final_nuc_ab[-1]

    one_cycle_cyt = np.array(mult_cycles_cyt[len(mult_cycles_cyt)-num_datapoints+num_nans:])
    one_cycle_nuc = np.array(mult_cycles_nuc[len(mult_cycles_cyt)-num_datapoints+num_nans:])

    # plotting
    plotting.plot_abundances(final_tspan[:-num_nans], one_cycle_cyt, one_cycle_nuc)
    plotting.plot_volume_ratio(t_range * time_scalar, nuc_vols, cell_vols)
    plotting.plot_abundance_ratio(final_tspan[:-num_nans], one_cycle_cyt, one_cycle_nuc)
    plotting.plot_concentration_ratio(final_tspan, one_cycle_cyt, one_cycle_nuc, cv_func, nv_func, num_nans)
    plotting.plot_multiple_cycles(final_tspan, mult_cycles_cyt, mult_cycles_nuc, num_cycles)


if __name__ == '__main__':
    main()
    sys.exit(0)

# NOTES
# From the volume analysis script: Average cycle duration: 14.270588235294118 frames, which is equal
# to 71.3529411764706 minutes. The data is interpolated to 100 datapoints, meaning that every datapoint is 0.716 minutes
