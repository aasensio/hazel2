import numpy as np
import matplotlib.pyplot as pl

tmp = np.loadtxt('rtcoef', skiprows=1)

fig, ax = pl.subplots(nrows=2, ncols=4, figsize=(15, 6), sharex=True)

N = 10.0**4.34
dz = 1000.0 * 1e5 / 10.0
for i in range(8):
    ax.flat[i].plot(tmp[:, 1], N * dz * tmp[:, i + 2])