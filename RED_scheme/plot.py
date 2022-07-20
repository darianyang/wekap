#!/usr/bin/env python
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, NullFormatter)
import numpy
import h5py

n_simulations = 1
scheme = '20-100conWE_lt36C2'
iterations = 500 

raw_datasets = [f'../{scheme}/raw_rates.npy']

datasets = [f'../{scheme}/rates.npy']

raw_rates = numpy.array(numpy.load(raw_datasets[0]))

rates = numpy.array(numpy.load(datasets[0]))

print(raw_rates, rates)

raw_rates_mean = numpy.mean(raw_rates, axis=0)
raw_rates_sem = numpy.std(raw_rates, axis=0)/numpy.sqrt(n_simulations)
rates_mean = numpy.mean(rates, axis=0)
rates_sem = numpy.std(rates, axis=0)/numpy.sqrt(n_simulations)

raw_rates_uci = raw_rates_mean + raw_rates_sem*1.96
raw_rates_lci = raw_rates_mean - raw_rates_sem*1.96
rates_uci = rates_mean + rates_sem*1.96
rates_lci = rates_mean - rates_sem*1.96

f, axarr = pyplot.subplots(1, sharex=True)
f.set_size_inches(8,6)
    
xs = numpy.arange(1,iterations+1,1)    

print(numpy.shape(xs), numpy.shape(raw_rates_mean))
#pyplot.plot(xs, rates)
#pyplot.show()

standard_rate, = axarr.semilogy(xs, raw_rates, c='dodgerblue', linewidth = 2, label="original WE")
#axarr.fill_between(xs, raw_rates_uci, raw_rates_lci, facecolor="dodgerblue", linewidth=0, alpha=0.4)
red_rate, = axarr.semilogy(xs, rates, c='red', linewidth=2, label="RED WE")
#axarr.fill_between(xs, rates_uci, rates_lci, facecolor="red", linewidth=0, alpha=0.4)

legend = axarr.legend(loc="lower right", handles=[red_rate,standard_rate],prop={'size':18})
legend.get_frame().set_linewidth(0)
for label in legend.get_lines():
    label.set_linewidth(3)
 
pyplot.xlabel('molecular time (ns)', size=22)
pyplot.ylabel('rate constant (s$^{-1}$)', fontsize=22)   

axarr.xaxis.set_major_locator(MultipleLocator(50))
axarr.xaxis.set_major_formatter(FormatStrFormatter('%d'))
axarr.xaxis.set_minor_locator(MultipleLocator(10))
axarr.tick_params(axis='x', which='major', length=10, width=0.75, labelsize=20)
axarr.tick_params(axis='x', which='minor', length=5, width=0.75)

locmaj = matplotlib.ticker.LogLocator(base=10.0, numticks=15)
axarr.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(2,4,6,8), numticks=15)
axarr.yaxis.set_minor_locator(locmin)
axarr.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axarr.tick_params(axis='y', which='major', length=5, width=0.75, labelsize=20)
axarr.tick_params(axis='y', which='minor', length=2.5, width=0.5)
for label in axarr.yaxis.get_ticklabels()[1::2]:
  label.set_visible(False)

#pyplot.xlim(0,200)
#axarr.set_ylim(bottom=1, top=1e6)
pyplot.tight_layout()
pyplot.savefig(f'../{scheme}/rates.pdf')

print("Original WE:", raw_rates_mean,"+/-", raw_rates_uci)
print("RED WE:", rates_mean, "+/-", rates_uci)
