import sys

import scipy.interpolate

import plotting
import numpy as np
import pandas as pd
from math import log
from scipy.integrate import odeint
from scipy.interpolate import interp1d

########################
### Model parameters ###
########################
kc = 0.25  # synthesis rate of protein in cytoplasm
kd = log(2) / 35  # degradation rate for protein
kIn = log(2) / 10  # rate of translocation into nucleus
kOut = log(2) / 10  # rate of translocation out of nucleus
peak_of_nuc_vol = 81  # position of the latest local high of the nuclear volume
nuc_div_tp = 91  # simulated point at which nuclear division takes place

num_cycles = 6  # number of cycles to include in the multiple-cycles plot
num_datapoints = 200  # the desired number of datapoints that is solved for within the time-axis

########################
### Global variables ###
########################
cell_vols, nuc_vols, nuc_surface_areas = [], [], []
a_func, cv_func, nv_func = interp1d.__class__, interp1d.__class__, interp1d.__class__

                                            ########################
                                            ### Data preparation ###
                                            ########################

def load_and_adjust_data():
    """
    Takes the by the script exported averages and alters these based on various assumptions:
      - Nuclear division takes place approximately 10 time-steps before the end (which represents whole-cell splitting)
      - Nuclear division duration is very short, and thus assumed to be instant (split off at one time-point)
      - Nuclear volume and surf area stays constant from local high up until division point
    :return:
    """
    global cell_vols, nuc_vols, nuc_surface_areas
    averages = pd.read_excel("./averages.xlsx")

    # nuclear volume over time (altered to simulate instant nuclear division)
    nv = averages.nuc_volumes.values
    nv[peak_of_nuc_vol:nuc_div_tp] = nv[peak_of_nuc_vol]  # up until nuc div, the nuc vol stays constant at its high
    nv[nuc_div_tp:] = [nv[0]] * len(nv[nuc_div_tp:])  # after nuclear separation, the volume returns to the starting val
    nuc_vols = np.append(nv, [nv[0], nv[0]])

    # nuclear surface area over time (altered such that it simulates a proportional decrease compared to nuclear volume)
    # basically the same as done above
    nsa = averages.nuc_surface_areas.values
    nsa[peak_of_nuc_vol:nuc_div_tp] = nsa[peak_of_nuc_vol]
    nsa[nuc_div_tp:] = [nsa[0]] * len(nsa[nuc_div_tp:])
    nuc_surface_areas = np.append(nsa, [nsa[0], nsa[0]])

    # whole-cell volumes
    cv = averages.cell_volumes.values
    cv[-1] = cv[0]  # we assume that the whole-cell volume returns to start value
    # TODO this isn't the case, the volume lost at division is only ~15%, not ~30%
    cell_vols = np.append(cv, [cv[0], cv[0]])

    # generate interp1d functions so the area and vols can be retrieved for any t
    define_interpolated_functions()


def define_interpolated_functions():
    """
    Creates scipy interp1d/interpolated functions based on the averaged nuclear surface, nuclear volume and
    whole-cell volume. Using these, the area and vols can be retrieved for any t within the 0 to 100 range.
    :return:
    """
    global a_func, cv_func, nv_func
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
                                                ### The Model ###
                                                #################

def dp_dt(y, t, k_d, k_s, k_in, k_out):
    """
    Model that attempts to describe the change in cellular and nuclear abundances of a protein using ODEs. It is a
    constitutive model, i.e. it provides a relationship between the behavior of a protein and the forces acting on it.
    No form of regulation is taken into account in this model.

    Equation description:
       dC/dt: (1) synthesis of protein (2) net transfer into the cyt (prop to nuc surf area) (3) deg. of prot in cyt
       dN/dt: (1) net transfer into the nucleus (proportional to nuc surf. area) (2) degradation of prot in nucleus

    :param y: array holding the previous predictions of the nuclear and cytoplasmic abundances
    :param t: the new timestep for which predictions are made
    :param k_d: degradation rate of protein
    :param k_s: synthesis rate of protein
    :param k_in: nuclear import rate of protein
    :param k_out: nuclear export rate of protein
    :return:
    """
    C, N = y  # cellular and nuclear abundance of protein
    A = a_func(t)  # nuclear surface area at t
    Vc = cv_func(t)  # whole cell volume at t
    Vn = nv_func(t)  # nuclear volume at t
    k_in = k_in / np.average(nuc_surface_areas)  # import rate scaling using the average nuclear surface area
    k_out = k_out / np.average(nuc_surface_areas)  # export rate scaling using the average nuclear surface area


    dC_dt = k_s * 20 + A * (-k_in * C / (Vc - Vn) + k_out * N / Vn) - k_d * C
    dN_dt = A * (k_in * C / (Vc - Vn) - k_out * N / Vn) - k_d * N

    return [dC_dt, dN_dt]


                                            ###########################
                                            ### Running simulations ###
                                            ###########################

def simulate(cp0, np0):
    """
    Function that simulates the protein abundance dynamics over the duration of a cell cycle based on the initial
    conditions and the desired amount of cycles.
    :param cp0: initial cytoplasmic protein abundance
    :param np0: initial nuclear protein abundance
    :return:
    """
    nuc_ab_loss_frac = (np.amax(nuc_vols) - nuc_vols[-1]) / np.amax(nuc_vols)
    whole_vol_loss_frac = 1 - cell_vols[-1] / cell_vols[-4]  # fraction of volume that is lost to the daughter bud
    print(f"Percentage of nuclear protein abundance lost at nuclear division: {round(nuc_ab_loss_frac * 100, 2)}%")
    print(f"Percentage of whole-cell volume lost at division (end of cycle): {round(whole_vol_loss_frac * 100, 2)}%")

    # arrays to store the results of multiple cycles in
    mult_cycles_cyt, mult_cycles_nuc = [], []

    # generate split t_span to be able to simulate nuclear division event at 'nuc_div_tp'
    tspan_whole = np.linspace(0, 99, num_datapoints)
    closest = find_nearest(tspan_whole, nuc_div_tp)
    tspan_before, tspan_after = np.split(tspan_whole, np.where(tspan_whole == closest)[0] + 1)

    # simulate certain amount of cycles
    for i in range(num_cycles):

        # solve the ode up until the nuclear division event
        sols_before = odeint(dp_dt, [cp0, np0], tspan_before, args=(kd, kc, kIn, kOut))
        sols_bf_ev_cyt = sols_before[:, 0]
        sols_bf_ev_nuc = sols_before[:, 1]

        # simulate nuclear division event (reduction of nuclear abundance proportional to loss in nuc vol)
        sols_after_event = np.array([
            sols_bf_ev_cyt[-1],  # nothing happens to the cellular abundance at nuclear split
            (1 - nuc_ab_loss_frac) * sols_bf_ev_nuc[-1]  # remove the calc. percentage from the ab.
        ])

        # use the altered initial conditions to simulate the dynamics after the nuclear division event
        sols_after = odeint(dp_dt, sols_after_event, tspan_after, args=(kd, kc, kIn, kOut))
        sols_after_cyt = sols_after[:, 0][~np.isnan(sols_after[:, 0])]
        sols_after_nuc = sols_after[:, 1][~np.isnan(sols_after[:, 1])]

        # simulate the bud separation event by reducing the cytoplasmic abundance in proportion to the volume loss
        # percentage (~ 28%) determined using the volume analysis pipeline
        sols_after_cyt[-1] = (1 - whole_vol_loss_frac) * sols_after_cyt[-1]

        final_tspan = np.concatenate((tspan_before, tspan_after))  # merge the two split t-spans
        # according to the volume analysis script, the average duration of one cycle is 71.35 minutes.
        # let's modify the time-axis using this knowledge (for purpose of plotting on a real time axis)

        final_cyt_ab = np.concatenate((sols_before[:, 0], sols_after_cyt))
        final_nuc_ab = np.concatenate((sols_before[:, 1], sols_after_nuc))
        mult_cycles_cyt.extend(final_cyt_ab)
        mult_cycles_nuc.extend(final_nuc_ab)

        # final values of this simulation become the next simulation's initial conditions
        cp0, np0 = final_cyt_ab[-1], final_nuc_ab[-1]

    return final_tspan, mult_cycles_cyt, mult_cycles_nuc


def main():
    load_and_adjust_data()
    define_interpolated_functions()

    cp0 = 200  # initial cytoplasmic protein abundance
    np0 = 10  # initial nuclear protein abundance

    # perform model simulations
    final_tspan, mult_cycles_cyt, mult_cycles_nuc = simulate(cp0, np0)

    # get one cycle for plotting purposes
    one_cycle_cyt = np.array(mult_cycles_cyt[len(mult_cycles_cyt)-num_datapoints:])
    one_cycle_nuc = np.array(mult_cycles_nuc[len(mult_cycles_cyt)-num_datapoints:])

    # plotting
    params = {"kd": round(kd, 5), "kIn": round(kIn, 4), "kOut": round(kOut, 4)}
    plotting.plot_abundances(final_tspan, one_cycle_cyt, one_cycle_nuc, params)
    plotting.plot_concentration_ratio(final_tspan, one_cycle_cyt, one_cycle_nuc, cv_func, nv_func, params)
    plotting.plot_multiple_cycles(final_tspan, mult_cycles_cyt, mult_cycles_nuc, num_cycles, params)


if __name__ == '__main__':
    main()
    sys.exit(0)

# NOTES
# From the volume analysis script: Average cycle duration: 14.270588235294118 frames, which is equal
# to 71.3529411764706 minutes. The data is interpolated to 100 datapoints, meaning that every datapoint is 0.716 minutes
