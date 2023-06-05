import sys
import plotting
import numpy as np
import pandas as pd
from math import log, pow
from scipy.integrate import odeint
from scipy.interpolate import interp1d

averages_file = "./averages.xlsx"

########################
                                            ### Model parameters ###
                                            ########################

kc = 0.25  # synthesis rate of protein in cytoplasm
kd = log(2) / 35  # degradation rate for protein
kIn = log(2) / 10  # rate of translocation into nucleus
kOut = log(2) / 10  # rate of translocation out of nucleus
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
    averages = pd.read_excel(averages_file)

    # nuclear volume over time (altered to simulate instant nuclear division)
    nv = averages.nuc_volumes.values
    peak_of_nuc_vol = np.argmax(nv)  # position of the volume peak of the nucleus
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
    """
    Finds the value in the array that is closest to 'value'
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

                                                #################
                                                ### The Model ###
                                                #################

def get_k_in(t, k_in):
    """
    Function that increases the nuclear import rate in the first third of the cell cycle. This
    is done by usage of a gaussian distribution.
    :param t:
    :param k_in:
    :return:
    """
    k_in_adj = k_in / np.average(nuc_surface_areas)

    # return 1 / sqrt(2 * pi * 70) * exp(-(t - 14)**2 / 2 * 0.03) + k_in_adj
    return 1 / sqrt(2 * pi * 16) * exp(-(t - 14)**2 / 2 * 0.03) + k_in_adj



def get_k_out(t, k_out):
    """
    Function that linearly increases the nuclear export rate from time point 85 in order
    to simulate transiently increasing leakiness of the nucleus when approaching karyokinesis.
    :param t:
    :param k_out:
    :return:
    """
    k_out_adj = k_out / np.average(nuc_surface_areas)

    return 0.1 * pow(1.1, t - 100) + k_out_adj


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
    k_out = get_k_out(t, k_out)
    k_in = get_k_in(t, k_in)  # import rate scaling using the average nuclear surface area
    # k_in = k_in / np.average(nuc_surface_areas)
    # k_out = k_out / np.average(nuc_surface_areas)

    ts.append(t)
    kins.append(k_in)
    kouts.append(k_out)

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
        sols_be = odeint(dp_dt, [cp0, np0], tspan_before, args=(kd, kc, kIn, kOut))

        # simulate nuclear division event (reduction of nuclear abundance proportional to loss in nuc vol)
        ab_after_event = np.array([
            sols_be[:, 0][-1],  # nothing happens to the cellular abundance at nuclear split
            (1 - nuc_ab_loss_frac) * sols_be[:, 1][-1]  # remove the calc. percentage from the ab.
        ])

        # simulate the dynamics after the nuclear division event
        sols_ae = odeint(dp_dt, ab_after_event, tspan_after, args=(kd, kc, kIn, kOut))

        # simulate the bud separation event by reducing the cytoplasmic abundance in proportion to the volume loss
        # percentage (~ 28%) determined using the volume analysis pipeline
        sols_ae[:, 0][-1] = (1 - whole_vol_loss_frac) * sols_ae[:, 0][-1]

        final_cyt_ab = np.concatenate((sols_be[:, 0], sols_ae[:, 0]))
        final_nuc_ab = np.concatenate((sols_be[:, 1], sols_ae[:, 1]))

        # add this cycle's data to the rest of the data
        mult_cycles_cyt.extend(final_cyt_ab)
        mult_cycles_nuc.extend(final_nuc_ab)

        # final values of this simulation become the next simulation's initial conditions
        cp0, np0 = final_cyt_ab[-1], final_nuc_ab[-1]

    return tspan_whole, mult_cycles_cyt, mult_cycles_nuc


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
    plotting.plot_volumes(final_tspan, cv_func, nv_func)
    plotting.plot_abundances(final_tspan, one_cycle_cyt, one_cycle_nuc, params)
    plotting.plot_concentration_ratio(final_tspan, one_cycle_cyt, one_cycle_nuc, cv_func, nv_func, params)
    plotting.plot_multiple_cycles(final_tspan, mult_cycles_cyt, mult_cycles_nuc, num_cycles, params)


if __name__ == '__main__':
    main()
    sys.exit(0)

# NOTES
# From the volume analysis script: Average cycle duration: 14.270588235294118 frames, which is equal
# to 71.3529411764706 minutes. The data is interpolated to 100 datapoints, meaning that every datapoint is 0.716 minutes
