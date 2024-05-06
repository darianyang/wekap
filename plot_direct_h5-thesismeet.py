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

from scipy.stats import hmean
from scipy.stats.mstats import gmean
#from scipy.signal import savgol_filter

plt.style.use("/Users/darian/github/wedap/wedap/styles/default.mplstyle")
#plt.style.use("/Users/darian/github/wedap/styles/poster.mplstyle")

# TODO: function to make 4 panel plot
    # Plot of P_A, P_B, rate_AB, rate_BA, all as function of WE iteration

# TODO: transition this to CLI program like mdap and wedap (eventually maybe make mkap)
class Kinetics:

    def __init__(self, h5="direct.h5", min_state=None, max_state=None, prefix=None, statepop="direct",
                 tau=100e-12, state=1, label=None, units="rates", ax=None, savefig=False, color=None):
        """ TODO
        Parameters
        ----------
        h5 : str
            Name of output file from WESTPA w_direct or w_ipa, e.g. `direct.h5`
        tau : float
            The resampling interval of the WE simualtion.
            This should be in seconds, default 100ps = 100 * 10^-12 (s).
        state : int
            State for flux calculations (flux into `state`), 0 = A and 1 = B.
        label : str
            Data label.
        statepop : str
            'direct' for state_population_evolution from direct.h5 or
            'assign' for labeled_populations from assign.h5.
        units : str
            Can be `rates` (default) or `mfpts`.
        ax : mpl axes object
        savefig : bool
            Wether or not to save the figure to a png.
        """
        # read in direct.h5 file
        self.direct_h5 = h5py.File(h5, "r")
        # temp solution for getting assign.h5, eventually can just get rid of it
        # since I don't think the color/labeled population is as useful
        self.assign_h5 = h5py.File(h5[:-9] + "assign.h5", "r")

        self.prefix = prefix
        self.tau = tau
        self.state = state
        self.label = label
        self.units = units
        self.statepop = statepop
        self.color = color

        self.min_state = min_state
        self.max_state = max_state
        if min_state and max_state:
            self.states = [i for i in range(min_state, max_state, 1)]
            self.n_states = len(self.states)
        else:
            self.states = None
            self.n_states = 1

        if self.statepop == "direct":
            # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
            self.state_pops = np.array(self.direct_h5["state_pop_evolution"])
            # state A = label 0, state B = label 1
            self.state_pop_a = np.array([expected[2] for expected in self.state_pops[:,0]])
            self.state_pop_b = np.array([expected[2] for expected in self.state_pops[:,1]])
        elif self.statepop == "assign":
            # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
            self.state_pops = np.array(self.assign_h5["labeled_populations"])
            # state A = label 0, state B = label 1
            #self.state_pop_a = np.array([expected[0] for expected in self.state_pops[:,0]])
            #self.state_pop_a = np.sum(self.state_pops[:,0], axis=1)
            #print(state_pop_a)
            #np.sum(self.state_pop_a)
            self.state_pop_a = np.sum(self.state_pops[:,0], axis=1)
            self.state_pop_b = np.sum(self.state_pops[:,1], axis=1)

        # create new fig or plot onto existing
        if ax is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax = ax
            self.fig = plt.gcf()

        self.savefig = savefig

    def extract_rate(self):
        """
        Get the raw rate array from one direct.h5 file.
        
        Returns
        -------
        rate_ab, ci_lb_ab, ci_ub_ab
        """

        # flux evolution dataset from cumulative evolution mode:
        # When calculating time evolution of rate estimates, 
        # ``cumulative`` evaluates rates over windows starting with --start-iter and 
        # getting progressively wider to --stop-iter by steps of --step-iter.
        fluxes = np.array(self.direct_h5["target_flux_evolution"])

        # conditional fluxes are macrostate to macrostate
        # 2 dimensions: [(0 -> 0, 0 -> 1), 
        #                (1 -> 0, 1 -> 1)] 
        # I want 0 -> 1
        #fluxes = np.array(h5["conditional_flux_evolution"])[:,:,1]

        # third column (expected) of the state (A(0) or B(1)) flux dataset (flux into state b = 1)
        flux_ab = np.array([expected[2] for expected in fluxes[:,self.state]])
        # CIs in rate (s^-1) format (divided by tau)
        ci_lb_ab = np.array([expected[3] for expected in fluxes[:,self.state]]) * (1/self.tau)
        ci_ub_ab = np.array([expected[4] for expected in fluxes[:,self.state]]) * (1/self.tau)

        # TODO: update to be a state assignment attr
        # norm by state pop A if calculating A --> B
        # if self.state == 1:
        #     state_pop = self.state_pop_a
        # # norm by state pop B if calculating B --> A
        # elif self.state == 0:
        #     state_pop = 1 - self.state_pop_a
        # # TODO: temp fix
        # else:
        #     state_pop = self.state_pop_a

        # assign the state of target flux flow
        if self.state == 1:
            state_pop = self.state_pop_a
        elif self.state == 0:
            state_pop = self.state_pop_b
        else:
            print("Currently only support state 0 or state 1.")

        # 2 different approaches here, can norm by state_pop_a (sum of weights in a)
        # but since 2 state system, could also use 1 - state_pop_b since all not in b are in a
        flux_ab = flux_ab / state_pop
        #flux_ab = flux_ab / state_pop_a
        #flux_ab = flux_ab / (1 - state_pop_b)

        # convert from tau^-1 to seconds^-1
        rate_ab = flux_ab * (1/self.tau)
        
        return rate_ab, ci_lb_ab, ci_ub_ab

    def plot_rate(self, title=None, moltime=True):
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
        if moltime:
            # multiply by tau (ps)
            iterations *= 100
            # convert to ns
            iterations = np.divide(iterations, 1000)

        if self.units == "mfpts":
            mfpt_ab = 1 / rate_ab
            self.ax.plot(iterations, mfpt_ab, label=self.label)
            #ax.fill_between(iterations, mfpt_ab - (1/ci_lb_ab), mfpt_ab + (1/ci_ub_ab), alpha=0.5)
            self.ax.set_ylabel("MFPT ($s$)")
        elif self.units == "rates":
            self.ax.plot(iterations, rate_ab, color=self.color)#, label=self.label)
            self.ax.fill_between(iterations, rate_ab - ci_lb_ab, rate_ab + ci_ub_ab, alpha=0.5,
                                 label=self.label, color=self.color)
            self.ax.set_ylabel("Rate Constant ($s^{-1}$)")

        if moltime:
            self.ax.set_xlabel(r"Molecular Time (ns)")
        else:
            self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
        
        self.ax.set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        self.ax.set_title(title)

        return rate_ab

    def plot_statepop(self):
        """
        Plot the state populations
        """
        # # divide k_AB by P_A for equilibrium rate correction (AB and BA steady states)
        # state_pops = np.array(self.direct_h5["state_pop_evolution"])
        # state A = label 0, state B = label 1
        #state_pop_a = np.array([expected[2] for expected in self.state_pops[:,0]])
        #state_pop_b = np.array([expected[2] for expected in self.state_pops[:,1]])
        # state_pop_a = np.sum(self.state_pops[:,0], axis=1)
        # state_pop_b = np.sum(self.state_pops[:,1], axis=1)

        # WE iterations
        iterations = np.arange(0, len(self.state_pop_a), 1)

        # plot both state population evolutions
        self.ax.plot(iterations, self.state_pop_a, label="State A")
        self.ax.plot(iterations, self.state_pop_b, label="State B")
        self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")
        self.ax.set_ylabel("State Population")

        #return self.state_pop_a, self.state_pop_b

    def plot_exp_vals(self, ax=None, f_range=False, d2d1=False, f_range_all=False):
        """
        f_range : bool
            Set to True to use mark 25-67 s^-1 as the k_D1D2 rate.
        d2d1 : bool
            Set to True to also include k_D2D1.
        """
        if ax is None:
            ax = self.ax
        if self.units == "rates":
            if f_range_all:
                # ax.axhline(60, alpha=1, color="tab:orange", label="4F k$_{D1D2}$", ls="--")
                # ax.axhline(25, alpha=1, color="tab:green", label="7F k$_{D1D2}$", ls="--")
                # ax.axhline(60, alpha=1, color="tab:red", ls="--")
                # ax.axhline(25, alpha=1, color="tab:green", ls="--")
                ax.axhline(60, alpha=1, color="tab:orange", ls="--")
                ax.axhline(25, alpha=1, color="tab:green", ls="--")
            elif f_range:
                # DTY 19F rates of 25-60 for k_D1D2
                ax.axhspan(25, 60, alpha=0.25, color="grey", label="NMR k$_{D1D2}$")
                if d2d1:
                    ax.axhspan(135, 179, alpha=0.25, color="tan", label="NMR k$_{D2D1}$")
            else:
                # D1-->D2 ~ 20-50, D2-->D1 ~ 100-150
                ax.axhline(25, color="k", ls="--", label="k$_{D2D1}$")
                if d2d1:
                    ax.axhline(150, color="red", ls="--", label="k$_{D1D2}$")
        elif self.units == "mfpts":
            # converted to mfpt = 1 / rate
            ax.axhline(1/150, color="k", ls="--", label="MFPT$_{D2D1}$")
            if d2d1:
                ax.axhline(1/25, color="red", ls="--", label="MFPT$_{D1D2}$")
        else:
            raise ValueError(f"You put {self.units} for unit, which must be `mfpts` or `rates`.") 

class KineticsMulti(Kinetics):
    # TODO: make this update scheme more automatically
    # TODO: prob best to convert to child class for multiple direct.h5 files
        # or is there a better way here? should think about this

    def __init__(self, h5s=["direct.h5"], labels=None, tau=10**-10, state=1, units="rates", 
                 ax=None, savefig=False):
        """
        For multiple direct.h5 files.

        Parameters
        ----------
        h5 : list or strs
            List of multiple direct.h5 files.
        labels : list of strs
            List of labels cooresponding to each direct.h5 file.
        """
        super().__init__(tau=tau, state=state, units=units, ax=ax, savefig=savefig)
        self.direct_h5s = h5s
        self.labels = labels
        if ax is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax = ax
            self.fig = plt.gcf()

    def plot_multi_def_rates(self, def_subplot=None):
        """
        Show individual rate plots for multiple state definitions.
        """
        final_rates = []
        for i, direct_h5 in enumerate(self.direct_h5s):
            self.h5 = direct_h5
            if self.labels:
                self.label = self.labels[i]
            rates = self.plot_rate()
            # add 2 item list: angle | final rate value
            final_rates.append(rates[-1])

        plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        #plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        # TODO: separate into different methods
        # plot the rates at various state definitions
        if def_subplot:
            ax2 = self.fig.add_subplot(def_subplot)
            ax2.plot(self.states, final_rates, color="k")
            ax2.set_xlabel("Angle State Definition", fontsize=11, labelpad=4, color="grey")
            self.plot_exp_vals(ax2)
            ax2.tick_params(colors="grey")
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        #self.ax.set_ylim(10**-9, 10**10)

    def plot_std_error_rate_reps(self, iterations=500, reps=3, def_subplot=None):
        """
        ### get rates with std errors for multiple state definitions
        Make a plot of multiple replicates and std err for each tstate def
            Iteration (X) VS Rate (Y) with std error
        And maybe later a plot of multiple state defs (X) vs rate(Y) with error
            Maybe both in one plot again?
        """
        rate_evolutions = np.zeros(shape=(iterations, self.n_states, reps))
        final_rate_scalers = np.zeros(shape=(self.n_states, reps))
        ogscheme = self.scheme
        # for v00, v01, v02, (threshv00)
        for vi, ver in enumerate(range(0,reps)):
        #for vi, ver in enumerate(["v00", "v01", "v02", "threshv00"]):
            for ai, angle in enumerate(self.states):
                self.scheme = f"{ogscheme}/{self.prefix}{angle}_v0{ver}"
                #print(self.scheme)
                self.label = f"> {angle}°"
                # calc and fill out rate
                rate_ab, ci_lb_ab, ci_ub_ab = self.extract_rate()
                # make 3d array of rates: rows = tau, columns = angle, depth = version
                rate_evolutions[:, ai, vi] = rate_ab[:iterations]
                #rate_evolutions[:, ai, vi] = rate_ab
                # array for final rate value
                final_rate_scalers[ai, vi] = rate_ab[-1]

        # plot avg and std error of rate evo
        avg_rate_evo = np.average(rate_evolutions, axis=2)
        stdev_rate_evo = np.std(rate_evolutions, axis=2)
        sterr_rate_evo = stdev_rate_evo / np.sqrt(rate_evolutions.shape[2])
        for ai, angle in enumerate(self.states):
            avg = avg_rate_evo[:,ai]
            err = sterr_rate_evo[:,ai]
            self.ax.plot(avg, label=f"< {angle}°")
            self.ax.fill_between(np.arange(0, iterations, 1), avg - err, avg + err, alpha=0.25)

        self.plot_exp_vals(self.ax)
        plt.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        # plt.legend(plot, [f"<{i}°" for i in states],
        #                 loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False)
        #plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
        plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        # state definition subplot
        if def_subplot:
            # plot the rates at various state definitions
            ax2 = self.fig.add_subplot(def_subplot)

            # plot avg and std error of final rates
            avg_final_scaler = np.average(final_rate_scalers, axis=1)
            stdev_final_scaler = np.std(final_rate_scalers, axis=1)
            #print(final_rate_scalers.shape[1])
            sterr_final_scaler = stdev_final_scaler / np.sqrt(final_rate_scalers.shape[1])

            # TODO: table of final avg +/- stderror values
            rate_table = np.zeros((len(self.states), 3))
            rate_table[:,0] = self.states
            rate_table[:,1] = avg_final_scaler
            rate_table[:,2] = sterr_final_scaler
            #np.savetxt("rate_table.csv", rate_table, delimiter="\t")
            import pandas as pd
            df = pd.DataFrame(rate_table, columns=["State", "Average", "Stderr"])
            print(df)

            ax2.plot(self.states, avg_final_scaler, color="k")
            ax2.fill_between(self.states, 
                             avg_final_scaler - sterr_final_scaler, 
                             avg_final_scaler + sterr_final_scaler, 
                             alpha=0.25, color="k")

            ax2.set_xlabel("Angle State Definition", fontsize=11, labelpad=4, color="grey")
            plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
            self.plot_exp_vals(ax2)
            ax2.tick_params(colors="grey")

        self.ax.set_ylabel("Rate Constant ($s^{-1}$)")
        self.ax.set_xlabel(r"WE Iteration ($\tau$=100ps)")

        return final_rate_scalers

    # To smooth the rate per tstate def plot
    def smooth(self, y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode="same")
        return y_smooth

    def plot_per_state_def(self, ax):
        final_rate_scalers = self.plot_std_error_rate_reps()

        # plot the rates at various state definitions        
        # plot avg and std error of final rates
        avg_final_scaler = np.average(final_rate_scalers, axis=1)
        stdev_final_scaler = np.std(final_rate_scalers, axis=1)
        sterr_final_scaler = stdev_final_scaler / np.sqrt(final_rate_scalers.shape[1])

        ax.plot(self.states, avg_final_scaler, color="k")
        ax.fill_between(self.states, 
                        avg_final_scaler - sterr_final_scaler, 
                        avg_final_scaler + sterr_final_scaler, 
                        alpha=0.25, color="k")

        #ax.set_xlabel("Angle State Definition", fontsize=16, labelpad=4)
        ax.set_xlabel("Angle State Definition")
        self.ax = ax
        self.plot_exp_vals(self.ax)
        #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])

        x = self.states
        y = avg_final_scaler
        # attempts at fitting a curve
        # window size 9, polynomial order 3
        #yhat = savgol_filter(y, 9, 3)
        #yhat = self.smooth(y, 10)

        # technically dx is constant so don't need diff(x)
        dydx = np.diff(y) / np.diff(x)
        ddydx = np.diff(dydx)
        print(ddydx)
        return dydx, ddydx

################################################
###### a more complicated derivative plot ######
################################################
#fig, ax = plt.subplots(figsize=(10,7), sharey=True, ncols=1)
#k = Kinetics(56,67,ax=ax)
#k.scheme = f"65tdef_1a43_c2_we/4b_64_v00"
#k.plot_rate()

#fig, ax = plt.subplots(figsize=(12,5), sharey=True, ncols=2)
#plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
#Kinetics(56,67,ax=ax).plot_multi_def_rates()
#Kinetics(56,66, ax=ax, state=1).plot_std_error_rate_reps(reps=3)
#Kinetics(26, 37, ax=ax).plot_multi_def_rates()
#dydx, ddydx = Kinetics(56,67,ax=ax[0]).plot_per_state_def(ax[1])

# fig2, ax2 = plt.subplots(ncols=1)
# angles = [i for i in range(56, 66, 1)]
# ax2.plot(angles[:-1], dydx, color="k")
# #ax2[1].plot(angles[:-2], ddydx)
# ax2.set_xlabel("Angle State Definition")
# ax2.set_ylabel("∆ Rate Constant ($s^{-1}$)")
# #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
# #ax2.set_yscale("symlog")
# ax2.set_ylim(-10000, 100)
# ax2.set_xlim(58, 65)
# fig2.tight_layout()

#plt.savefig("figures/2kod_c2_20-80_v00.png", dpi=300, transparent=True)

###################################
###### single direct.h5 plot ######
###################################
# k = Kinetics(scheme="o_and_c2_angle_2kod_v00")
# k.plot_rate()
# k.plot_exp_vals()

######################################
###### common multi state plots ######
######################################
#k = Kinetics(4, 13, scheme="2kod_oa_lt15oa", prefix="1d_", state=1, savefig=True)
#k.plot_multi_def_rates(ver="v00")
# k = Kinetics(45, 53, scheme="2kod_oa_65c2", prefix="2d_", state=1, savefig=True)
# k.plot_std_error_rate_reps(reps=3, def_subplot=326, title="2KOD >65° C2 Angle")

# for c2 in range(65, 66):
#     k = Kinetics(45, 53, scheme=f"2kod_oa_{c2}c2", prefix="2d_", state=1, savefig=False)
#     plt.ylim(10**-8, 10**7)
#     rates = k.plot_std_error_rate_reps(reps=3, def_subplot=326, title=f"2KOD >{c2}° C2 Angle")

#################################################
### save the 50° oa value and plot them later ###
#################################################
# c2_rates = []
# for c2 in range(60, 72):
#     k = Kinetics(50, 51, scheme=f"2kod_oa_{c2}c2", prefix="2d_", state=1, savefig=False)
#     rates = k.plot_std_error_rate_reps(reps=3, title=f"2KOD >{c2}° C2 Angle")
#     cols = [c2, np.average(rates), (np.std(rates)/np.sqrt(3))]
#     c2_rates.append(cols)

# c2_rates = np.array(c2_rates)
# print(c2_rates)
# fig3, ax3 = plt.subplots()

# print(c2_rates[:,0] - c2_rates[:,2])

# ax3.plot(c2_rates[:,0], c2_rates[:,1], color="k")
# ax3.fill_between(c2_rates[:,0], 
#                  c2_rates[:,1] - c2_rates[:,2], 
#                  c2_rates[:,1] + c2_rates[:,2], 
#                  alpha=0.25, color="k")

# ax3.set_xlabel("Angle State Definition")
# ax3.set_ylabel("Rate Constant ($s^{-1}$)")
# #plt.yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
# #plt.yscale("symlog", subs=[2, 3, 4, 5, 6, 7, 8, 9])
# plt.ylim(-10, 250)
# k.plot_exp_vals(ax3)
# #ax3.tick_params(colors="grey")
# plt.legend()
# plt.tight_layout()
# fig3.savefig("figures/c2_60-71.png", dpi=300, transparent=True)

#plt.tight_layout()
#plt.show()


# TODO: a heatmap plot with X definition, Y definition, and Z as the final rate constant
#########################
### state def heatmap ###
#########################
# TODO: make a custom cmap with the center around the experimental rate
# c2_rates = []
# for c2 in range(60, 72):
#     k = Kinetics(45, 54, scheme=f"2kod_oa_{c2}c2", prefix="2d_", state=1, savefig=False)
#     rates = k.plot_std_error_rate_reps(reps=3, title=f"2KOD >{c2}° C2 Angle")
#     #cols = [c2, np.average(rates), (np.std(rates)/np.sqrt(3))]
#     #print(rates[-1:].shape)
#     c2_rates.append(rates)

# x = [i for i in range(60, 72)]
# y = [i for i in range(45, 54)]
# z = np.average(np.array(c2_rates), axis=2)

# import matplotlib.colors

# fig4, ax4 = plt.subplots()
# mesh = ax4.pcolormesh(y, x, z, shading="auto", norm=matplotlib.colors.LogNorm(), cmap="afmhot")
# cbar = fig4.colorbar(mesh)
# cbar.set_label("Rate Constant ($s^{-1}$)")
# ax4.set(xlabel="Orientation Angle (°)", ylabel="C2 Angle (°)")
# fig4.tight_layout()
# fig4.savefig("figures/angle_heatmap.png", dpi=300)


######################
### 1d oapdt plots ###
######################
# fig, ax = plt.subplots()
# # systems = ["2kod_v00", "lo_pH_v00", "150end_v00"]
# #systems = ["2kod_v00"]
# systems = ["2kod_v00", "4F_v00", "7F_v00"]
# #names = ["2kod WT", "Low pH", "150-end"]
# #names = ["From WE"]
# for i in range(len(systems)):
#     #k = Kinetics(scheme=f"gt2500oapdt_gt70c2/{systems[i]}", state=1, label=names[i], ax=ax)
#     k = Kinetics(scheme=f"oapdt_c2_2dgrid/{systems[i]}/2500oapdt_70c2", state=0, label=systems[i], ax=ax)
#     k.plot_rate()
#     # k = Kinetics(scheme=f"gt2500oapdt_gt65c2/{systems[i]}", state=1, label=names[i], ax=ax)
#     # k.plot_rate()
# k.plot_exp_vals(f_range=True)
# plt.legend()
# plt.tight_layout()
# #plt.ylim(0, 1000)
# plt.show()
#fig.savefig("figures/wt_only_2500oapdt_70c2.png", dpi=300, transparent=True)

############################
### oapdt c2 plots/grids ###
############################
# TODO: the reverse rates look pretty good, explore this more
# fig, ax = plt.subplots(ncols=2, sharey=True, figsize=(10,5))
# #systems = ["2kod_v00", "lo_pH_v00", "150end_v00"]
# systems = ["2kod_v00"]
# #systems = ["lo_pH_v00"]
# final_rates = []
# for c2 in range(66,73):
# #for c2 in range(70,71):
#     #for oapdt in range(2100, 2700, 100):
#     for oapdt in range(2500, 2600, 100):
#         k = Kinetics(scheme=f"oapdt_c2_2dgrid/{systems[0]}/{oapdt}oapdt_{c2}c2",
#                     state=0, label=f"{str(c2)}_{oapdt}", ax=ax[0])
#         rate = k.plot_rate()
#         final_rates.append(rate[-1])
# ax[1].plot([i for i in range(66,73)], final_rates)
# #ax[1].plot([i for i in range(2100,2700,100)], final_rates)
# k.plot_exp_vals(f_range=True, ax=ax[0])
# #plt.legend()
# ax[0].legend()
# ax[1].set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
# #ax[1].set_xlabel("Orientation Angle Product")
# #ax[1].set_xlabel("C2 Angle (°)")
# plt.tight_layout()
# plt.show()
# fig.savefig("figures/wt_lo_150end_2500oapdt_70c2.png", dpi=300, transparent=True)
#fig.savefig("figures/testing_2500oapdt.png", dpi=300, transparent=True)
#fig.savefig("figures/testing_70c2.png", dpi=300, transparent=True)

# TODO: try/except block for OS error for failed kinetics 


##################################################
### oapdt c2 plots of multiple def final rates ###
##################################################
# fig, ax = plt.subplots(ncols=2, sharey=True, figsize=(10,5))
# #systems = ["2kod_v00", "lo_pH_v00", "150end_v00"]
# systems = ["2kod_v00"]
# #systems = ["lo_pH_v00"]
# final_rates = []
# for c2 in range(66,73):
# #for c2 in range(70,71):
#     #for oapdt in range(2100, 2700, 100):
#     for oapdt in range(2500, 2600, 100):
#         k = Kinetics(scheme=f"oapdt_c2_2dgrid/{systems[0]}/{oapdt}oapdt_{c2}c2")
#         rate = k.plot_rate()
#         final_rates.append(rate[-1])
# ax[1].plot([i for i in range(66,73)], final_rates)
# #ax[1].plot([i for i in range(2100,2700,100)], final_rates)
# k.plot_exp_vals(f_range=True, ax=ax[0])
# #plt.legend()
# ax[0].legend()
# ax[1].set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
# #ax[1].set_xlabel("Orientation Angle Product")
# ax[1].set_xlabel("C2 Angle (°)")
# plt.tight_layout()
# plt.show()
# fig.savefig("figures/wt_lo_150end_2500oapdt_70c2.png", dpi=300, transparent=True)

#######################
### oamax c2 plots  ###
#######################
# TODO: run some ssWE, need a good tstate
# fig, ax = plt.subplots()
# k = Kinetics(scheme=f"oamax_c2_2dgrid_2/WT_v00/51oamax_76c2", state=1, statepop="direct", ax=ax)
# k.plot_rate()
# # k = Kinetics(scheme=f"oamax_c2_2dgrid_2/WT_v00/51oamax_76c2", state=0, statepop="direct", ax=ax)
# # k.plot_rate()
# #k2 = Kinetics(scheme=f"oamax_c2_2dgrid/WT_v00/52oamax_76c2", state=0, statepop="direct", ax=ax)
# #k2.plot_rate()
# k.plot_exp_vals(f_range=True, ax=ax)
# #plt.legend()
# plt.tight_layout()
# plt.show()

#######################
### oamax c2 plots  ###
#######################
def hmean_stderr(multi_k):
    """
    From: The Standard Errors of the Geometric and Harmonic Means 
          and Their Application to Index Numbers
          Author(s): Nilan Norris
          Source: The Annals of Mathematical Statistics , 
                  Dec., 1940, Vol. 11, No. 4 (Dec., 1940), pp. 445-448

    s_H = (1/alpha^2)(s_(1/x)/sqrt(n-1))
    alpha = 1 / H
    s_(1/x) = stdev(1/x)
                  
    Parameters
    ----------
    multi_k : 2darray
        Array of n arrays.

    Returns
    -------
    hmeans : 1darray
        Harmonic means.
    s_H : 1darray
        Standard error of harmonic means.
    """
    #print(len(multi_k))
    hmeans = hmean(multi_k, axis=0)
    #print(hmeans.shape)
    alpha = np.reciprocal(hmeans)
    #print(alpha.shape)
    std_recip = np.std(np.reciprocal(multi_k), axis=0)
    #print(std_recip.shape)
    s_H = (1 / (alpha**2)) * (std_recip / np.sqrt(len(multi_k) - 1))
    #print(s_H.shape)

    return hmeans, s_H

fig, ax = plt.subplots()
# k = Kinetics(scheme=f"oa1_oa2_c2/WT_v00/60oamax_72c2", state=3, statepop="direct", ax=ax)
# k.plot_rate()
# k = Kinetics("1a43_oamin_oa12/3d_v00/45oa/direct.h5", state=1, statepop="direct", ax=ax)
# k.plot_rate()
# k.plot_exp_vals()
# k = Kinetics("1a43_oamin_oa12/3d_v01/45oa/direct.h5", state=1, statepop="direct", ax=ax)
# k.plot_rate()

# TODO: these could be the args for the CLI/GUI, along with arg for error preference
# k = Kinetics("2kod_4dmin_v00/direct.h5", state=1, statepop="assign", ax=ax)
# k.plot_rate()
#k.plot_exp_vals()
#k.plot_statepop()

def plot_multi_run(sys="WT", ax=None, oa="56", color=None):
    """
    Plot multiple runs of a system.
    Use bayesian bootstrapping for error estimates.
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = plt.gcf()

    # multi rate list per replicate
    multi_k = []

    # calc and append the 1st item which is the rates
    for i in range(5):
        print(f"D1D2_lt16oa/{sys}_v0{i}/{oa}oa/direct.h5")
        k = Kinetics(f"D1D2_lt16oa/{sys}_v0{i}/{oa}oa/direct.h5", state=1, ax=ax)
        #k = Kinetics(f"oa1_oa2_c2/{sys}_v0{i}/{oa}oa_72c2/direct.h5", state=1, ax=ax)
        #k = Kinetics(f"oa1_oa2/{sys}_v0{i}/{oa}oa/direct.h5", state=1, ax=ax)
        multi_k.append(k.extract_rate()[0])
    # harmonic mean
    #multi_k_avg, multi_k_stderr = hmean_stderr(multi_k)
    # hmean and error from bayesian bootstrapping
    #multi_k_avg = hmean(multi_k)
    from error.bootstrap import bayboot_multi
    multi_k_stderr = bayboot_multi(multi_k, repeat=1000)
    print(multi_k_stderr.shape)

    # arithmetic mean
    multi_k_avg = np.average(multi_k, axis=0)
    #multi_k_std = np.std(multi_k, axis=0)
    #multi_k_stderr = multi_k_std / np.sqrt(len(multi_k))
    #multi_k_stderr = np.rot90(np.vstack((multi_k_stderr,multi_k_stderr)))
    #print(multi_k_stderr.shape)

    iterations = np.arange(0, len(multi_k_avg), 1)
    # multiply by tau (ps)
    iterations *= 100
    # convert to ns
    iterations = np.divide(iterations, 1000)
    ax.plot(iterations, multi_k_avg, color=color)
    # plt.fill_between(iterations, multi_k_avg - multi_k_stderr, 
    #                  multi_k_avg + multi_k_stderr, alpha=0.5, label="WT")
    ax.fill_between(iterations, multi_k_stderr[:,0],
                    multi_k_stderr[:,1], alpha=0.5, label=sys, color=color)
    #print(np.average(multi_k, axis=0).shape)
    #print("WT AVG and STDERR: ", multi_k_avg[-1], multi_k_stderr[-1])

    #print(multi_k_avg[-1])
    return k, multi_k_avg[-1], multi_k_stderr[-1], multi_k

# comparing means
# multi_k_avg = gmean(multi_k, axis=0)
# multi_k_stderr = bayboot_multi(multi_k, gmean, repeat=1000)
# plt.plot(iterations, multi_k_avg)
# plt.fill_between(iterations, multi_k_avg - multi_k_stderr[:,0], 
#                  multi_k_avg + multi_k_stderr[:,1], alpha=0.5, label="GMEAN")
# multi_k_avg = hmean(multi_k, axis=0)
# multi_k_stderr = bayboot_multi(multi_k, hmean, repeat=1000)
# plt.plot(iterations, multi_k_avg)
# plt.fill_between(iterations, multi_k_avg - multi_k_stderr[:,0], 
#                  multi_k_avg + multi_k_stderr[:,1], alpha=0.5, label="HMEAN")

# 19F
# k = Kinetics("oa1_oa2_c2/4F_v00/55oa_72c2/direct.h5", state=1, statepop="direct", ax=ax, label="4F", color="tab:red")
# print("4F AVG: ", k.plot_rate()[-1])
# CIs
# print(k.extract_rate()[1][-1])
#print(k.extract_rate()[2][-1])
# k = Kinetics("oa1_oa2_c2/7F_v00/55oa_72c2/direct.h5", state=1, statepop="direct", ax=ax, label="7F", color="tab:green")
# print("7F AVG: ", k.plot_rate()[-1])
#print(k.extract_rate()[2][-1])

# multi class test
# k = KineticsMulti([f"oa1_oa2_c2/WT_v0{i}/55oa_72c2/direct.h5" for i in range(5)], state=1, ax=ax)
# k.plot_multi_def_rates()
# finals = []
# for i in range(54, 59):
#     k, f = plot_multi_run("WT", ax=ax, oa=i)
#     finals.append(f)


angle = 12
#k, f = plot_multi_run("WT", ax=ax, oa=angle)
k, f1, e1, mk1 = plot_multi_run("WT", ax=ax, oa=angle, color="tab:blue")
k, f2, e2, mk2 = plot_multi_run("4F", ax=ax, oa=angle, color="tab:orange")
#k, f3, e3, mk3 = plot_multi_run("7F", ax=ax, oa=angle, color="tab:green")
# k.plot_exp_vals(f_range_all=True, ax=ax)
# plt.legend(loc=4, ncol=1)
# ax.set_yscale("log", subs=[2, 3, 4, 5, 6, 7, 8, 9])
# #plt.ylim(0,5000)
# plt.xlabel("Molecular Time (ns)")
# plt.ylabel("Rate Constant (s$^{-1}$)")
errors = np.flip(np.rot90(np.vstack((e1, [0,0], [0,0], e2, [3.34, 3.34], [8.26,8.26], [0,0], [0,0])), k=-1), axis=1)
#print(np.rot90(errors))
print(errors)
ax.axhline(60, alpha=1, color="tab:orange", ls="--")
ax.axhline(25, alpha=1, color="tab:green", ls="--")
ax.set_yscale("log")

#fig, ax = plt.subplots(figsize=(10,5))
fig, ax = plt.subplots(figsize=(12,5))
ax.bar(["WT-eqWE", "WT-ssWE", "MSM", "4F-eqWE", "4F-NMR", "4F-NMRr", "WT-ssWEr", "MSMr"], 
       [f1,         1.5,        0.015,    f2,      60,       134,       173,        0.05], 
       color=["tab:blue", "tab:blue", "grey", "tab:orange", "magenta", "magenta", "tab:blue", "grey"],
       yerr=errors, capsize=10)
ax.set_ylabel("Rate Constant (s$^{-1}$)")

def plot_mk_scatter(mk, ax, label="WT-WE"):
    mk = [i[-1] for i in mk]
    print(mk)
    ax.scatter([label for _ in range(5)], mk, color="k")
plot_mk_scatter(mk1, ax, "WT-eqWE")
plot_mk_scatter(mk2, ax, "4F-eqWE")
#plot_mk_scatter(mk3, ax, "7F-WE")


ax.set_yscale("log")
#plt.xticks(fontweight="bold")
#plt.xticks(rotation=30, ha='right')
plt.tight_layout()
#plt.show()
#plt.savefig("figures/wt_mean_comp.png", dpi=300, transparent=True)
#plt.savefig("figures/bar_all_poster.png", dpi=300, transparent=True)
plt.savefig("figures/05May2024_all_attempts_so_far.png", dpi=300, transparent=True)
