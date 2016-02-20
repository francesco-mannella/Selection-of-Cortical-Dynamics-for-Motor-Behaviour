import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

data = np.loadtxt("activation_data")
w = np.loadtxt("weights")

fig = plt.figure(figsize=(12,4))
gs = gridspec.GridSpec(4,12)
ax1 = fig.add_subplot(gs[:2,:8])
ax2 = fig.add_subplot(gs[2:,:8])
ax3 = fig.add_subplot(gs[:,8:])

ax1.plot(data)

ax2.imshow(data.T, aspect="auto", interpolation="none")
ax2.set_axis_off()

e = np.linalg.eigvals(w)
rho = np.max(np.abs(e))
ax3.scatter(np.real(e), np.imag(e) )
ax3.set_xlim([-rho, rho])
ax3.set_ylim([-rho, rho])
ax3.set_xlabel("real")
ax3.set_ylabel("imag")
plt.tight_layout()
plt.show()

