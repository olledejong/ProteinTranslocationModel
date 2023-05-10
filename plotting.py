import matplotlib.pyplot as plt
import numpy as np
plt.style.use('seaborn-v0_8-dark')

time_scalar = 0.7135  # scalar based on average duration of cycle to scale back to real minute axis

def plot_abundances(tspan, y1, y2, params):
    fig, ax1 = plt.subplots()
    fig.suptitle(f"Cytoplasmic and nuclear protein abundances over time\nParams: {params}")
    ax1.set_xlabel('Time (minutes)')
    ax1.grid(False)
    ax1.set_ylabel("Cytoplasmic protein abundance", color='orange')
    ax1.plot(tspan * time_scalar, y1, color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='darkred')
    ax2.plot(tspan * time_scalar, y2, color='darkred')
    ax2.tick_params(axis='y', labelcolor='darkred')
    plt.show()


def plot_concentration_ratio(final_tspan, one_cycle_cyt, one_cycle_nuc, cv_func, nv_func, params):
    # get cell / nuc volumes for all timepoints
    cell_vols = [cv_func(t).flatten()[0] for t in final_tspan]
    nuc_vols = [nv_func(t).flatten()[0] for t in final_tspan]
    cyt_vols = [i - j for i, j in zip(cell_vols, nuc_vols)]

    # with the abundances, the event is distinctly happening over one time-step, this is not the case for the volumes.
    # this is because the volumes are first manipulated to represent nuclear division at one time-step, but after that
    # they are interpolated using the interp1d method. This causes a more transient simulation of that division.
    # therefore we need to slightly alter these to prevent plotting artifacts.
    nuc_vols[181:184] = [nuc_vols[180]] * len(nuc_vols[181:184])  # TODO make this dynamic (non-hardcoded indexes)
    cyt_vols[-2] = cyt_vols[-3]

    # calculate the concentrations
    c_con = [i / j for i, j in zip(one_cycle_cyt.tolist(), cyt_vols)]
    n_con = [i / j for i, j in zip(one_cycle_nuc.tolist(), nuc_vols)]

    con_ratio = [i / j for i, j in zip(n_con, c_con)]

    plt.plot(final_tspan * time_scalar, c_con, c='darkred', lw=2)
    plt.title(f"Cytoplasmic protein concentration\nParams: {params}")
    plt.xlabel("Time (minutes)")
    plt.ylabel("Concentration")
    plt.show()

    plt.plot(final_tspan * time_scalar, n_con, c='darkred', lw=2)
    plt.title(f"Nuclear protein concentration\nParams: {params}")
    plt.xlabel("Time (minutes)")
    plt.ylabel("Concentration")
    plt.show()

    plt.plot(final_tspan * time_scalar, con_ratio, c='darkred', lw=2)
    plt.title(f"Nuclear to cytoplasmic protein concentration ratio\nParams: {params}")
    plt.xlabel("Time (minutes)")
    plt.ylabel("Ratio")
    plt.show()


def plot_multiple_cycles(final_tspan, cyt_ab_cycles, nuc_ab_cycles, num_cycles, params):
    t_axis = np.linspace(0, final_tspan[-1] * num_cycles, len(cyt_ab_cycles))
    fig, ax1 = plt.subplots()
    fig.suptitle(f"Cytoplasmic and nuclear protein abundances over time (multiple cycles)\nParams: {params}")
    ax1.set_xlabel('Time (minutes)')
    ax1.grid(False)
    ax1.set_ylabel("Cytoplasmic protein abundance", color='orange')
    ax1.plot(t_axis * time_scalar, cyt_ab_cycles, color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='darkred')
    ax2.plot(t_axis * time_scalar, nuc_ab_cycles, color='darkred')
    ax2.tick_params(axis='y', labelcolor='darkred')
    plt.show()
