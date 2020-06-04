import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# plt.switch_backend("Agg")

seeds = np.array([1612835610, 1809, 1947713535, 2226841396, 2801740628, 315775938, 3213595428, 3454752513, 4157420972, 4260127168, 724315229])

ind = np.array([4, 13, 36])

kparal = np.load("data/default/power_spectrum.npz", allow_pickle=True)["kparal"][0]
kparal = kparal[kparal>0]

some_ps = np.zeros((len(kparal), 11))

indx = ind[2]

for ii in range(len(seeds)):

	random_seed = seeds[ii]
	filename = "data/%i/power_spectrum.npz" %random_seed

	lc_data = np.load(filename, allow_pickle=True)

	P = lc_data["power_spectrum"]
	kperp = lc_data["kperp"]
	# kparal = lc_data["kparal"][0]

	P = P[1:,int(len(kparal))+2:]
	# kparal = kparal[kparal>0]
	kperp = kperp[kperp>0]

	some_ps[:,ii] = P[indx,:]


fig, ax = plt.subplots(1, 1, sharey=True, gridspec_kw={"hspace":0.15, "wspace":0.05, "left":0.15})#, subplot_kw={"xscale":"log", "yscale":'log'})
for ii in range(len(seeds)):
	ax.plot(kparal, (some_ps[:,ii] - np.mean(some_ps, axis=1))/np.mean(some_ps, axis=1))
	# ax[1].plot(kparal, np.log10(P[ind[1],:]))
	# ax[2].plot(kparal, np.log10(P[ind[2],:]))

ax.set_title(r"$k_{\perp}=$%.3f [Mpc$^{-1} h$]" %kperp[indx], size= 10)
# ax[1].set_title(r"$k_{\perp}=$%.3f [Mpc$^{-1} h$]" %kperp[ind[1]], size = 10)
# ax[2].set_title(r"$k_{\perp}=$%.3f [Mpc$^{-1} h$]" %kperp[ind[2]], size=10)

ax.set_ylabel(r"$[P(k)-\bar{P}(k)] / \bar{P}(k)$", fontsize =10)
ax.set_xlabel(r"$k_{\parallel}$ [Mpc$^{-1} h$]", fontsize =12)

fig.savefig("figures/variability_seeds-newlimit-3", dpi=150)

# plt.figure()
# plt.imshow(np.log10(P).T, vmin=2, vmax=10, origin='lower', aspect='auto', extent=(np.min(kperp[kperp!=0]), np.max(kperp), np.min(kparal[kparal[0]>0]), np.max(kparal[kparal[0]>0])))
# plt.colorbar(label=r"mK$^2$Mpc$^3$h$^{-3}$")
# plt.xlabel("1/Mpc")
# plt.ylabel("1/Mpc")
# plt.xscale("log")
# plt.yscale("log")
# plt.savefig("figures/power_spectrum2d-newlimit")
# plt.close()