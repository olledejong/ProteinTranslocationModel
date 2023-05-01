import sys
import plotting
import numpy as np
import pandas as pd
from math import log
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline, interp1d, UnivariateSpline

peak_of_nuc_vol = 81  # the latest local high of the nuclear volume
nuc_div_tp = 91  # simulated point at which nuclear division takes place


#########################
### Data preparations ###
#########################

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def load_data():
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

    return cv, nv, nsa


###############################
### The model & simulations ###
###############################


def dp_dt(y, t, kd, kc, k_in, k_out):
    c, n = y  # cellular and nuclear abundance of protein

    # equations that describe the change in cellular and nuclear abundances
    # dC/dt: (1) synthesis of protein (2) net transfer into the cyt (prop to nuc surf area) (3) deg. of prot in cyt
    # dN/dt: (1) net transfer into the nucleus (proportional to nuc surf. area) (2) degradation of prot in nucleus
    dC_dt = kc * 50 + a_func(t) * (-k_in * c / (cv_func(t) - nv_func(t)) + k_out * n / nv_func(t)) - kd * c
    dN_dt = a_func(t) * (k_in * c / (cv_func(t) - nv_func(t)) - k_out * n / nv_func(t)) - kd * n

    return [dC_dt, dN_dt]


cell_vols, nuc_vols, nuc_surface_areas = load_data()
# interpolated vol/area functions for time dependent retrieval of values
t_range = np.arange(len(nuc_surface_areas))
a_func = interp1d(t_range, nuc_surface_areas, kind='linear', bounds_error=False)
cv_func = interp1d(t_range, cell_vols, kind='linear', bounds_error=False)
nv_func = interp1d(t_range, nuc_vols, kind='linear', bounds_error=False)


def main():
    # log(2)/10 is the rate of translocation, this is scaled by dividing it by the average nuclear surface
    kc = 0.6  # synthesis rate of protein in cytoplasm
    kd = log(2)/35  # degradation rate for protein
    kIn = log(2)/10 / np.average(nuc_surface_areas)  # rate of translocation into nucleus
    kOut = log(2)/10 / np.average(nuc_surface_areas)  # rate of translocation out of nucleus

    # initial conditions
    cp0 = 1400
    np0 = 55

    # generate split t_span to be able to simulate nuclear division event at 'nuc_div_tp'
    tspan_whole = np.linspace(0, 99, 200)
    closest = find_nearest(tspan_whole, nuc_div_tp)
    tspan_before, tspan_after = np.split(tspan_whole, np.where(tspan_whole == closest)[0] + 1)

    # solve the ode up until the event
    sols_before = odeint(dp_dt, [cp0, np0], tspan_before, args=(kd, kc, kIn, kOut))

    # simulate event (reduce nuclear abundance)
    perc_to_remove = (np.amax(nuc_vols) - nuc_vols[-1]) / np.amax(nuc_vols)  # reduction proportional to loss in nuc vol
    sols_after_event = np.array([
        sols_before[:, 0][-1],  # nothing happens to the cellular abundance at nuclear split
        sols_before[:, 1][-1] - (sols_before[:, 1][-1] * perc_to_remove)  # remove the calc. percentage from the ab.
    ])

    # simulate part after the event
    sols_after = odeint(dp_dt, sols_after_event, tspan_after, args=(kd, kc, kIn, kOut))
    sols_after_cyt_ab = sols_after[:, 0][~np.isnan(sols_after[:, 0])]
    sols_after_nuc_ab = sols_after[:, 1][~np.isnan(sols_after[:, 1])]

    # simulate the bud separation event by reducing the cytoplasmic abundance in proportion to the volume loss
    sols_after_cyt_ab[-1] = (1 - 0.28125901660197855) * sols_after_cyt_ab[-1]

    # some nans are produced because the outcome of the n_func / cv_func / nv_func is unknown at the last few timepoints
    num_nans = np.count_nonzero(np.isnan(sols_after[:, 0]))

    # final values (concatenated)
    final_tspan = np.concatenate((tspan_before, tspan_after))
    final_cyt_ab = np.concatenate((sols_before[:, 0], sols_after_cyt_ab))
    final_nuc_ab = np.concatenate((sols_before[:, 1], sols_after_nuc_ab))

    # plotting
    plotting.plot_abundances(final_tspan[:-num_nans], final_cyt_ab, final_nuc_ab)
    plotting.plot_volume_ratio(t_range, nuc_vols, cell_vols)
    plotting.plot_abundance_ratio(final_tspan[:-num_nans], final_cyt_ab, final_nuc_ab)
    plotting.plot_concentration_ratio(final_tspan, final_cyt_ab, final_nuc_ab, cv_func, nv_func, num_nans)


if __name__ == '__main__':
    main()
    sys.exit(0)
