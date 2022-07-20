"""
Plot the fluxes and rates from direct.h5 files.

list(h5):
['arrivals', 'avg_color_probs', 'avg_conditional_fluxes', 'avg_rates', 
'avg_state_probs', 'avg_total_fluxes', 'color_prob_evolution', 'conditional_arrivals',
'conditional_flux_evolution', 'conditional_fluxes', 'duration_count', 'durations', 
'rate_evolution', 'state_labels', 'state_pop_evolution', 'target_flux_evolution', 'total_fluxes']

  /target_flux_evolution [window,state]
    Total flux into a given macro state based on
    windows of iterations of varying width, as in /rate_evolution.
  /conditional_flux_evolution [window,state,state]
    State-to-state fluxes based on windows of
    varying width, as in /rate_evolution.

The structure of these datasets is as follows:
  iter_start
    (Integer) Iteration at which the averaging window begins (inclusive).
  iter_stop
    (Integer) Iteration at which the averaging window ends (exclusive).
  expected
    (Floating-point) Expected (mean) value of the observable as evaluated within
    this window, in units of inverse tau.
  ci_lbound
    (Floating-point) Lower bound of the confidence interval of the observable
    within this window, in units of inverse tau.
  ci_ubound
    (Floating-point) Upper bound of the confidence interval of the observable
    within this window, in units of inverse tau.
  stderr
    (Floating-point) The standard error of the mean of the observable
    within this window, in units of inverse tau.
  corr_len
    (Integer) Correlation length of the observable within this window, in units
    of tau.
"""

from cProfile import label
from typing import final
import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.style.use("/Users/darian/github/wedap/wedap/styles/default.mplstyle")

def plot_rate(scheme, tau=10**-10, state=1, label=None, units="rate", ax=None):
    """
    Plot the rate constant = target flux evolution AB / P_A 

    Parameters
    ----------
    scheme : str
        Scheme/directory for the direct.h5 file.
    tau : float
        The resampling interval of the WE simualtion.
        This should be in seconds, default 100ps = 100 * 10^-12 = 10^-10
    state : int
        State for flux calculations, 0 = A and 1 = B.
    label : str
        Data label.
    units : str
        Can be `rate` (default) or `mfpt`.
    ax : mpl axes object

    Returns
    -------
    rate_ab : ndarray
        Array of rates from A -> B in seconds^-1.
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = plt.gcf()

    # read in direct.h5 file
    h5 = h5py.File(f"{scheme}/direct.h5", "r")

    # flux evolution dataset from cumulative evolution mode:
    # When calculating time evolution of rate estimates, 
    # ``cumulative`` evaluates rates over windows starting with --start-iter and 
    # getting progressively wider to --stop-iter by steps of --step-iter.
    fluxes = np.array(h5["target_flux_evolution"])

    # conditional fluxes are macrostate to macrostate
    # 2 dimensions: [(0 -> 0, 0 -> 1), 
    #                (1 -> 0, 1 -> 1)] 
    # I want 0 -> 1
    #fluxes = np.array(h5["conditional_flux_evolution"])[:,:,1]

    # third column (expected) of the state (A/0 or B/1) flux dataset (flux into state b = 1)
    flux_ab = np.array([expected[2] for expected in fluxes[:,state]])
    # CIs
    ci_lb_ab = np.array([expected[3] for expected in fluxes[:,state]]) * (1/tau)
    ci_ub_ab = np.array([expected[4] for expected in fluxes[:,state]]) * (1/tau)
    # WE iterations
    iterations = np.arange(0, len(flux_ab), 1)

    # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
    state_pops = np.array(h5["state_pop_evolution"])
    # state A = label 0, state B = label 1
    state_pop_a = np.array([expected[2] for expected in state_pops[:,0]])
    #state_pop_b = np.array([expected[2] for expected in state_pops[:,1]])
    # 2 different approaches here, can norm by state_pop_a (sum of weights in a)
    # but since 2 state system, could also use 1 - state_pop_b since all not in b are in a
    flux_ab = flux_ab / state_pop_a
    #flux_ab = flux_ab / (1 - state_pop_b)

    # convert from tau^-1 to seconds^-1
    rate_ab = flux_ab * (1/tau)

    if units == "mfpt":
        mfpt_ab = 1 / rate_ab
        ax.plot(iterations, mfpt_ab, label=label)
        #ax.fill_between(iterations, mfpt_ab - (1/ci_lb_ab), mfpt_ab + (1/ci_ub_ab), alpha=0.5)
        ax.set_ylabel("MFPT ($s$)")
    elif units == "rate":
        ax.plot(iterations, rate_ab, label=label)
        ax.fill_between(iterations, rate_ab - ci_lb_ab, rate_ab + ci_ub_ab, alpha=0.5)
        ax.set_ylabel("Rate Constant ($s^{-1}$)")

    ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
    plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

    return rate_ab

def plot_statepop(scheme, ax=None):
    """
    Plot the state populations

    Parameters
    ----------
    scheme : str
        Scheme/directory for the direct.h5 file.
    ax : mpl axes object
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = plt.gcf()

    # read in direct.h5 file
    h5 = h5py.File(f"{scheme}/direct.h5", "r")

    # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
    state_pops = np.array(h5["state_pop_evolution"])
    # state A = label 0, state B = label 1
    state_pop_a = np.array([expected[2] for expected in state_pops[:,0]])
    state_pop_b = np.array([expected[2] for expected in state_pops[:,1]])

    # WE iterations
    iterations = np.arange(0, len(state_pop_a), 1)

    # plot both state population evolutions
    ax.plot(iterations, state_pop_a, label="State A")
    ax.plot(iterations, state_pop_b, label="State B")
    ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
    ax.set_ylabel("State Population")

    #return state_pop_a, state_pop_b

# TODO: function to make 4 panel plot
# Plot of P_A, P_B, rate_AB, rate_BA, all as function of WE iteration

fig, ax = plt.subplots(figsize=(10,6), sharey=True)

state = 1
units = "rate"

#plot_rate("conWEx_c2_rmsbbx", label="1A43 2D C2 RMSX", ax=ax, state=state)
# plot_rate("conWEx_c2_rmsbbx_750i", label="1A43 2D C2 RMSX", ax=ax, state=state)
# plot_rate("conWEx_rms_3.4", label="> 3.4A", ax=ax, state=state)
# plot_rate("conWEx_rms_3.5", label="> 3.5A", ax=ax, state=state)
# plot_rate("conWEx_rms_3.6", label="> 3.6A", ax=ax, state=state)
# plot_rate("20-100conWE_lt32C2", label="2KOD < 32° C2", ax=ax, state=state, units=units)

# plot_statepop("conWEx_c2_58", ax=ax)
# plot_statepop("conWEx_c2_60", ax=ax)

### get rates for multiple state definitions (WE c2x 5b)
# final_rates = []
# for angle in range(50, 64, 2):
#     rates = plot_rate(f"conWEx_c2_{angle}", label=f"> {angle}°", ax=ax, state=state, units=units)
#     # add 2 item list: angle | final rate value
#     final_rates.append([angle, rates[-1]])

### get rates for multiple state definitions (WE c2x 4b)
final_rates = []
for angle in range(58, 64, 2):
    rates = plot_rate(f"conWEx_4b/conWEx_c2_4b_{angle}", label=f"> {angle}°", ax=ax, state=state, units=units)
    # add 2 item list: angle | final rate value
    final_rates.append([angle, rates[-1]])

#plot_statepop("conWEx_c2_rmsbbx", ax=ax)
#plot_statepop("20-100conWE_lt32C2", ax=ax)

#plot_rate("20-100conWE_lt30C2", label="< 30°", ax=ax, state=state, units=units)
#plot_rate("20-100conWE_lt32C2", label="< 32°", ax=ax, state=state, units=units)
# plot_rate("20-100conWE_lt34C2", label="< 34°", ax=ax, state=state, units=units)
# plot_rate("20-100conWE_lt36C2", label="< 36°", ax=ax, state=state, units=units)

# plot_rate("20-100conWE_to65C2", label="to 65°", ax=ax, state=state, units=units)
# plot_rate("20-100conWE_to80C2", label="to 80°", ax=ax, state=state, units=units)
# plot_rate("20-100conWE_to82C2", label="to 82°", ax=ax, state=state, units=units)

#plot_rate("20-100conWE_lt31C2", label="< 31°", ax=ax, state=state)

# plot_rate("multi2kod_lt32C2", label="multi-2kod", ax=ax, state=state, units=units)
# plot_rate("v02_2kod_lt32C2", label="v02-2kod", ax=ax, state=state, units=units)

# shaded line from 100-150 s^-1 for experimental rate
# ax.axhline(100, color="k", ls="--")
# ax.axhline(150, color="k", ls="--", label="exp")

# D1-->D2 ~ 20-50, D2-->D1 ~ 100-150
ax.axhline(150, color="k", ls="--", label="k$_{D2D1}$")
ax.axhline(25, color="red", ls="--", label="k$_{D1D2}$")

# converted to mfpt = 1 / rate
# ax.axhline(1/100, color="k", ls="--")
# ax.axhline(1/150, color="k", ls="--", label="exp")

plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
#plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
#plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
fig.tight_layout()
#plt.savefig("WE_con_20-100_rates.png", dpi=300, transparent=True)

# plot the rates at various state definitions
# final_rates = np.array(final_rates)
# ax2 = fig.add_subplot(426)
# ax2.plot(final_rates[:,0], final_rates[:,1], color="k")
# ax2.set_xlabel("Angle State Definition", fontsize=11, labelpad=4)
# plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

plt.show()