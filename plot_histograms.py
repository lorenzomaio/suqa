import numpy as np
import matplotlib.pyplot as plt
from sys import argv

if len(argv)<1:
    print(f"usage: python3 {argv[0]} <fname> <nqe>")
    exit(0)

meas=np.loadtxt(argv[1])
nqe=int(argv[2])
Emin=1.9                                #float(argv[3])
Emax=11.3                               #float(argv[4])

fig,ax=plt.subplots()

phys=np.loadtxt("D4_g0.8_true_lams_deg.dat")
ct, bins,_ = ax.hist(phys,bins=200,density=True,label="exact spectrum")
ax.hist(np.linspace(Emin,Emax,2**nqe),bins=200,label="QPE grid")
ax.hist(meas[:,0],bins=bins,label="measures")

ax.set_xlabel("$E$")
plt.legend()
plt.show()
