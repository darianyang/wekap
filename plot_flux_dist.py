import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.style.use("/Users/darian/github/wedap/wedap/styles/default.mplstyle")

succ = np.loadtxt("/Users/darian/github/wekap/eg5/mmab_wt_v01/succ.txt", skiprows=1)

#plt.bar([i for i in range(0, 480)], succ[:,2])

plt.scatter([i for i in range(0, 480)], succ[:,2], c=succ[:,3], s=10)
cbar = plt.colorbar()
cbar.set_label("ADP RMSD ($\AA$)")

#plt.hist(succ[:,2]*1e100, 100)

#import seaborn as sns
#sns.jointplot(x=[i for i in range(0, 480)], y=succ[:,2])#, color=succ[:,3])

plt.yscale("log")
#plt.xscale("log")

plt.ylabel("Weight")
plt.xlabel("Recycling Event")

plt.tight_layout()
#plt.show()
plt.savefig("figures/eg5_event_weight_dist.png", dpi=300, transparent=True)

