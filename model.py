import sys
import traceback
import numpy as np
import pandas as pd
from math import exp, log, log10
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import parameters as params
import plotting as plot

ts = []
kouts = []
kins = []

averages_file = "./averages.xlsx"
ref_trace_file = "Sfp1_WT.xlsx"
ref_trace_file_path = f"./reference_traces/{ref_trace_file}"

# -------------------- Globals

nsa_func, wc_v_func, cyt_v_func, nuc_v_func = None, None, None, None

# -------------------- Data preparation

def load_and_adjust_data():
    """
    Takes the by the script exported averages and alters these based on various assumptions:
      - Nuclear division takes place approximately 10 time-steps before the end (which represents whole-cell splitting)
      - Nuclear division duration is very short, and thus assumed to be instant (split off at one time-point)
      - Nuclear volume and surf area stays constant from local high up until division point
    :return:
    """
    averages = pd.read_excel(averages_file)

    # variable for all datatypes
    nuc_vols = averages.nuc_volumes.values
    nuc_surf_areas = averages.nuc_surface_areas.values
    cell_vols = averages.cell_volumes.values
    cyt_vols = np.subtract(averages.cell_volumes.values, averages.nuc_volumes.values)

    t_range = np.arange(len(averages.cell_volumes.values))

    # plot the raw average volume values
    plot.plot_volumes(t_range, cyt_vols, averages.nuc_volumes.values, "before")

    # nuclear volume over time (altered to simulate instant nuclear division)
    peak_of_nuc_vol = np.argmax(nuc_vols)  # position of the volume peak of the nucleus
    nuc_vols[peak_of_nuc_vol:params.nuc_div_tp] = nuc_vols[peak_of_nuc_vol]  # nuc vol stays constant until nuc division
    nuc_vols[params.nuc_div_tp:] = [nuc_vols[0]] * len(nuc_vols[params.nuc_div_tp:])  # after div, the volume is reset

    # nuclear surface area over time (altered such that it simulates a proportional decrease compared to nuclear volume)
    # basically the same as done above
    nuc_surf_areas[peak_of_nuc_vol:params.nuc_div_tp] = nuc_surf_areas[peak_of_nuc_vol]
    nuc_surf_areas[params.nuc_div_tp:] = [nuc_surf_areas[0]] * len(nuc_surf_areas[params.nuc_div_tp:])

    # we assume that the whole-cell and cytosolic volume return to start value
    cell_vols[-1] = cell_vols[0]
    cyt_vols[-1] = cyt_vols[0]

    # plot the volumes after the alterations
    plot.plot_volumes(t_range, cyt_vols, nuc_vols, "after")


    return cell_vols, cyt_vols, nuc_vols, nuc_surf_areas


def define_interpolated_functions(cell_vols, cyt_vols, nuc_vols, nuc_surf_areas):
    """
    Creates scipy interp1d/interpolated functions based on the averaged nuclear surface, nuclear volume and
    whole-cell volume. Using these, the area and vols can be retrieved for any t within the 0 to 100 range.
    :return:
    """
    global nsa_func, wc_v_func, cyt_v_func, nuc_v_func
    # interpolated vol/area functions for time dependent retrieval of values
    t_range = np.arange(len(nuc_surf_areas))
    nsa_func = interp1d(t_range, nuc_surf_areas, kind='linear', bounds_error=False)
    wc_v_func = interp1d(t_range, cell_vols, kind='linear', bounds_error=False)
    cyt_v_func = interp1d(t_range, cyt_vols, kind='linear', bounds_error=False)
    nuc_v_func = interp1d(t_range, nuc_vols, kind='linear', bounds_error=False)

# -------------------- Other functions

def calc_concentration_ratio(final_tspan, one_cycle_cyt, one_cycle_nuc):
    """
    Function that calculates the nuclear-to-cytosolic concentration ratio of the protein
    :param final_tspan:
    :param one_cycle_cyt:
    :param one_cycle_nuc:
    :return:
    """
    # get cell / nuc volumes for all timepoints
    cytosolic_vols = np.array([cyt_v_func(t).flatten()[0] for t in final_tspan])
    nuclear_vols = np.array([nuc_v_func(t).flatten()[0] for t in final_tspan])

    # calculate the concentrations
    cyt_con = np.divide(one_cycle_cyt, cytosolic_vols)
    nuc_con = np.divide(one_cycle_nuc, nuclear_vols)

    return cyt_con, nuc_con, np.divide(nuc_con, cyt_con)

# -------------------- The model

def get_k_in(t, k_in, k_in_mp):
    """
    Function that can be used to define a nuclear import rate curve. In the report, a single Gaussian-like increase
    during G1, but also an increase during G1 ánd a decrease just before mitosis have been used for this purpose.
    However, this function might take any form.
    :param t:
    :param k_in:
    :param k_in_mp:
    :return:
    """
    return k_in_mp / 80 * exp(-(t - 12)**2 / 2 * 0.09) + k_in


def get_k_out(t, k_out, k_out_mp):
    """
    Function that can be used to define a nuclear export rate curve. In the report, both an exponential and a Gaussian
    function have been used for this purpose, but this might take any form.
    :param t:
    :param k_out:
    :param k_out_mp:
    :return:
    """
    return k_out_mp / 2 * exp(0.28 * (t - 115)) + k_out

def dp_dt(y, t, k_d, k_s, k_in, k_in_mp, k_out, k_out_mp, average_nuclear_surface_area):
    """
    Model that attempts to describe the change in cellular and nuclear abundances of a protein using ODEs. It is a
    constitutive model, i.e. it provides a relationship between the behavior of a protein and the forces acting on it.
    No form of regulation is taken into account in this model.

    Equation description:
       dC/dt: (1) synthesis of protein (2) net transfer into the cyt (prop to nuc surf area) (3) deg. of prot in cyt
       dN/dt: (1) net transfer into the nucleus (proportional to nuc surf. area) (2) degradation of prot in nucleus

    :param y: array holding the previous predictions of the nuclear and cytosolic abundances
    :param t: the new timestep for which predictions are made
    :param k_d: degradation rate of protein
    :param k_s: synthesis rate of protein
    :param k_in: nuclear import rate of protein
    :param k_out: nuclear export rate of protein
    :param k_in_mp: nuclear import rate multiplier
    :param k_out_mp: nuclear export rate multiplier
    :param average_nuclear_surface_area: the average nuclear surface area over the cell cycle
    :return:
    """
    if t > 99.0: t = 99.0  # prevent that t overshoots value of t (results in NaN values)

    C, N = y  # cellular and nuclear abundance of protein
    A = nsa_func(t)  # nuclear surface area at t
    Vc = wc_v_func(t)  # whole cell volume at t
    Vn = nuc_v_func(t)  # nuclear volume at t

    # scale using average nuclear surface
    k_in = k_in / average_nuclear_surface_area
    k_out = k_out / average_nuclear_surface_area

    k_out = get_k_out(t, k_out, k_out_mp)
    k_in = get_k_in(t, k_in, k_in_mp)

    ts.append(t)
    kins.append(k_in)
    kouts.append(k_out)

    dC_dt = k_s * 20 + A * (-k_in * C / (Vc - Vn) + k_out * N / Vn) - k_d * C
    dN_dt = A * (k_in * C / (Vc - Vn) - k_out * N / Vn) - k_d * N

    return [dC_dt, dN_dt]

# -------------------- Model simulation

def simulate(cell_vols, nuc_vols, nuc_surf_areas, kIn_base, kIn_mp, kOut_base, kOut_mp):
    """
    Function that simulates the protein abundance dynamics over the duration of a cell cycle based on the initial
    conditions and the desired amount of cycles.
    :return:
    """
    cp0, np0 = params.cp0, params.np0
    average_nuc_surf_area = np.average(nuc_surf_areas)

    nuc_ab_loss_frac = (np.amax(nuc_vols) - nuc_vols[-1]) / np.amax(nuc_vols)  # nuc ab loss prop. to nuc volume los
    whole_vol_loss_frac = 1 - cell_vols[-1] / cell_vols[-4]  # cytosolic ab loss proportional to whole-cell vol loss

    # arrays to store the results of multiple cycles in
    mult_cycles_cyt, mult_cycles_nuc = [], []

    # generate split t_span to be able to simulate nuclear division event at 'nuc_div_tp'
    tspan_whole = np.linspace(0, 99, params.num_datapoints)
    tspan_before, tspan_after = np.split(tspan_whole, np.where(tspan_whole == params.nuc_div_tp)[0])

    # simulate certain amount of cycles
    for i in range(params.num_cycles):

        # solve the ode up until the nuclear division event
        sols_be = odeint(dp_dt, [cp0, np0], tspan_before, args=(
            params.kd, params.ks, kIn_base, kIn_mp, kOut_base, kOut_mp, average_nuc_surf_area
        ))

        # simulate nuclear division event (reduction of nuclear abundance proportional to loss in nuc vol)
        ab_after_event = np.array([
            sols_be[:, 0][-1],  # nothing happens to the cellular abundance at nuclear split
            (1 - nuc_ab_loss_frac) * sols_be[:, 1][-1]  # remove the calc. percentage from the ab.
        ])

        # simulate the dynamics after the nuclear division event
        sols_ae = odeint(dp_dt, ab_after_event, tspan_after, args=(
            params.kd, params.ks, kIn_base, kIn_mp, kOut_base, kOut_mp, average_nuc_surf_area
        ))

        # simulate the bud separation event by reducing the cytosolic abundance in proportion to the volume loss
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

# -------------------- Main (start-point)

def main():
    cell_vols, cyt_vols, nuc_vols, nuc_surf_areas = load_and_adjust_data()
    define_interpolated_functions(cell_vols, cyt_vols, nuc_vols, nuc_surf_areas)

    ref_trace_pd = pd.read_excel(ref_trace_file_path)
    ref_trace = np.insert(ref_trace_pd.values.flatten(), 0, ref_trace_pd.columns[0])

    prev_mse = 1000  # ignore the assigned value. is reassigned every time a prediction by the model has a lower mse
    best_result = 'placeholder'  # will be filled with parameters for simulation that is most comparable to reference
    best_plot = 'placeholder'
    num_samples = 5  # number of samples to grab from logistic range of base values

    kin_base_samples = [x for x in 10 ** np.linspace(log10(log(2)), log10(log(2) / 4), num_samples)]
    kout_base_samples = [x for x in 10 ** np.linspace(log10(log(2) / 5), log10(log(2) / 20), num_samples)]

    count = 1

    for kout_base in kout_base_samples:
        for kin_base in kin_base_samples:
            for kin_to_multiply in kin_base_samples:
                for kout_to_multiply in kout_base_samples:
                    multiplier_kin = (kin_base_samples.index(kin_to_multiply) + 1)
                    multiplier_kout = (kout_base_samples.index(kout_to_multiply) + 1)

                    # perform model simulations
                    final_tspan, mult_cycles_cyt, mult_cycles_nuc = simulate(
                        cell_vols, nuc_vols, nuc_surf_areas,
                        kin_base, multiplier_kin, kout_base, multiplier_kout
                    )

                    # get final cycle for plotting purposes
                    one_cycle_cyt = np.array(mult_cycles_cyt[len(mult_cycles_cyt)-params.num_datapoints:])
                    one_cycle_nuc = np.array(mult_cycles_nuc[len(mult_cycles_cyt)-params.num_datapoints:])

                    # plotting
                    plot.plot_rates(ts, kouts, kins, count)
                    ts.clear()
                    kouts.clear()
                    kins.clear()

                    cyt_con, nuc_con, con_ratio = calc_concentration_ratio(final_tspan, one_cycle_cyt, one_cycle_nuc)

                    mse = plot.plot_prediction_vs_reference(
                        final_tspan, con_ratio, count, kin_base, kout_base, multiplier_kin, multiplier_kout, ref_trace
                    )
                    if mse < prev_mse:
                        best_result = {
                            'kIn base': kin_base,
                            'kOut base': kout_base,
                            'kIn mp': multiplier_kin,
                            'kOut mp': multiplier_kout
                        }
                        prev_mse = mse
                        best_plot = count
                    count += 1

    print(f"Done. Result with lowest Mean Squared Error ({prev_mse}) is: \n{str(best_result)}.\nSee plot {best_plot}.")


if __name__ == '__main__':
    try:
        main()
        sys.exit(0)
    except Exception:
        print("\nSomething went wrong while running the model:\n", traceback.format_exc())
        sys.exit(1)
