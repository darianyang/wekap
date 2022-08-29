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

import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.style.use("/Users/darian/github/wedap/wedap/styles/default.mplstyle")

# TODO: gotta convert to a class object

# TODO: function to make 4 panel plot
    # Plot of P_A, P_B, rate_AB, rate_BA, all as function of WE iteration

class Kinetics:

    def __init__(self, scheme=None, tau=10**-10, state=1, label=None, units="rates", ax=None):
        """
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
            Can be `rates` (default) or `mfpts`.
        ax : mpl axes object
        """
        self.scheme = scheme
        self.tau = tau
        self.state = state
        self.label = label
        self.units = units

        # create new fig or plot onto existing
        if ax is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax = ax
            self.fig = plt.gcf()

    def extract_rate(self):
        """
        Get the raw rate array from one direct.h5 file.
        """
        # read in direct.h5 file
        h5 = h5py.File(f"{self.scheme}/direct.h5", "r")

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
        flux_ab = np.array([expected[2] for expected in fluxes[:,self.state]])
        # CIs in rate (s^-1) format (divided by tau)
        ci_lb_ab = np.array([expected[3] for expected in fluxes[:,self.state]]) * (1/self.tau)
        ci_ub_ab = np.array([expected[4] for expected in fluxes[:,self.state]]) * (1/self.tau)

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
        rate_ab = flux_ab * (1/self.tau)
        
        return rate_ab, ci_lb_ab, ci_ub_ab

    def plot_rate(self, title=None):
        """
        Plot the rate constant = target flux evolution AB / P_A 

        Returns
        -------
        rate_ab : ndarray
            Array of rates from A -> B in seconds^-1.
        """
        rate_ab, ci_lb_ab, ci_ub_ab = self.extract_rate()

        # WE iterations
        iterations = np.arange(0, len(rate_ab), 1)

        if self.units == "mfpts":
            mfpt_ab = 1 / rate_ab
            self.ax.plot(iterations, mfpt_ab, label=self.label)
            #ax.fill_between(iterations, mfpt_ab - (1/ci_lb_ab), mfpt_ab + (1/ci_ub_ab), alpha=0.5)
            self.ax.set_ylabel("MFPT ($s$)")
        elif self.units == "rates":
            self.ax.plot(iterations, rate_ab, label=self.label)
            self.ax.fill_between(iterations, rate_ab - ci_lb_ab, rate_ab + ci_ub_ab, alpha=0.5)
            self.ax.set_ylabel("Rate Constant ($s^{-1}$)")

        self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        self.ax.set_title(title)

        return rate_ab

    def plot_statepop(self):
        """
        Plot the state populations
        """
        # read in direct.h5 file
        h5 = h5py.File(f"{self.scheme}/direct.h5", "r")

        # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
        state_pops = np.array(h5["state_pop_evolution"])
        # state A = label 0, state B = label 1
        state_pop_a = np.array([expected[2] for expected in state_pops[:,0]])
        state_pop_b = np.array([expected[2] for expected in state_pops[:,1]])

        # WE iterations
        iterations = np.arange(0, len(state_pop_a), 1)

        # plot both state population evolutions
        self.ax.plot(iterations, state_pop_a, label="State A")
        self.ax.plot(iterations, state_pop_b, label="State B")
        self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
        self.ax.set_ylabel("State Population")

        return state_pop_a, state_pop_b

    def plot_exp_vals(self):
        if self.units == "rates":
            # D1-->D2 ~ 20-50, D2-->D1 ~ 100-150
            self.ax.axhline(150, color="k", ls="--", label="k$_{D2D1}$")
            self.ax.axhline(25, color="red", ls="--", label="k$_{D1D2}$")
        elif self.units == "mfpts":
            # converted to mfpt = 1 / rate
            self.ax.axhline(1/150, color="k", ls="--", label="MFPT$_{D2D1}$")
            self.ax.axhline(1/25, color="red", ls="--", label="MFPT$_{D1D2}$")
        else:
            raise ValueError(f"You put {self.units} for unit, which must be `mfpts` or `rates`.") 

    def plot_multi_def_rates(self, ver="v00"):
        ### get rates for multiple state definitions (WE c2x 4b)
        final_rates = []
        for angle in range(56, 68, 1):
            self.scheme = f"1a43_c2_we/4b_{angle}_{ver}"
            self.label = f"> {angle}°"
            rates = self.plot_rate(title=f"1A43 {ver}")
            # add 2 item list: angle | final rate value
            final_rates.append([angle, rates[-1]])

        self.plot_exp_vals()

        plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        #plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        # TODO: separate into different methods
        # plot the rates at various state definitions
        final_rates = np.array(final_rates)
        ax2 = self.fig.add_subplot(528)
        ax2.plot(final_rates[:,0], final_rates[:,1], color="k")
        ax2.set_xlabel("Angle State Definition", fontsize=11, labelpad=4)
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        self.fig.tight_layout()
        #plt.savefig(f"figures/1a43_we_20-100_rates_{ver}.png", dpi=300, transparent=True)
        

    def plot_std_error_rate_reps(self):
        """ TODO (might be easier to update to class structure first)
        Make a plot of multiple replicates and std err for each tstate def
            Iteration (X) VS Rate (Y) with std error
        And maybe later a plot of multiple state defs (X) vs rate(Y) with error
            Maybe both in one plot again?
        """
        ### get rates for multiple state definitions (WE c2x 4b)
        final_rates = []
        for angle in range(56, 68, 1):
            rates = self.plot_rate(f"1a43_c2_we/4b_{angle}", label=f"> {angle}°", ax=self.ax, state=self.state, units=self.units, title=f"1A43")
            # add 2 item list: angle | final rate value
            final_rates.append([angle, rates[-1]])

        plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        #plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        self.fig.tight_layout()
        #plt.savefig("WE_con_20-100_rates.png", dpi=300, transparent=True)

        # plot the rates at various state definitions
        final_rates = np.array(final_rates)

        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        #plt.savefig(f"figures/1a43_we_20-100_rates.png", dpi=300, transparent=True)
        plt.show()


Kinetics().plot_multi_def_rates()
plt.show()