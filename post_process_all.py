import h5py
import numpy as np
import astropy.constants as const
from astropy.cosmology import Planck15 as cosmo
import astropy.units as un
from conversion import Lnu_BH_radio, Lnu_BH_quasar, Lnu_SFR, M_Eddington, fluxDensity, scaledFluxDensity, mpctoRadian, redshifttoFrequency, mKtoJy_per_sr
import matplotlib.cm as cm
from powerbox import get_power
import time
import os
import sys
import scipy.interpolate
from scipy import signal
from read_meraxes import ReadHeader, ReadGalaxiesSnap, ReadGridData

def indicesRandomization(position_dimension, bins_pos, distance):

	max_val = np.max(bins_pos) - distance

	start_pos = np.random.choice(position_dimension[position_dimension<=max_val])

	end_pos = start_pos + distance

	return (position_dimension>=start_pos) & (position_dimension<=end_pos)

minfrequency = 150e6
maxfrequency = 160e6
boxsize = 500/cosmo.h*un.Mpc
Ngrid = 256
startsnap = 60
endsnap = 61
max_Jy = 10e-3
bins_pos = np.linspace(0, boxsize.value, Ngrid+1)

random_seed = int(sys.argv[1])
inputdatadir = "/fred/oz009/bnasirud/meraxes/output/output-Genesis_L500_N2160-%i/" %random_seed
outputdir = "/fred/oz009/bnasirud/process/data/"
lcgridfilename = "meraxes_grids_48.hdf5"

# os.makedirs("data/%i" %random_seed,exist_ok=True)

galdesiredfields = ["dt","Pos","BlackHoleMass","BlackHoleAccretedHotMass","BlackHoleAccretedColdMass","Spin","Sfr"]
lcdesiredfields = ["LightconeBox","lightcone-z"]

hdffile = h5py.File(inputdatadir + lcgridfilename, "r")
redshift = hdffile["lightcone-z"][:]
hdffile.close()

# change redshift to frequency and only save some of the lightcone within 10 MHz band
frequency_lc_all = redshifttoFrequency(redshift)
lc_slice_sel=(frequency_lc<=maxfrequency) & (frequency_lc>=minfrequency)
redshift_lc = redshift_lc[lc_slice_sel]
frequency_lc = frequency_lc[lc_slice_sel]

nfreq = len(frequency_lc)
lightcone_fg = np.zeros((Ngrid,Ngrid,nfreq)) # set it to be big over z direction

#Read the number of cores and redshift from the main meraxes file
ncores,redshift,totngalaxies = ReadHeader(inputdatadir,startsnap,endsnap)


for isnap in range(startsnap,endsnap): 

	lightcone_fg_buff = lightcone_fg.tobytes()

	##########################################################################
	######################   GenerateLightconeFG  ###########################
	##########################################################################

	#Get the slice from the buffer so it can be populated
	# lightcone_fg = np.frombuffer(lightcone_fg_buff).reshape(Ngrid,Ngrid,len(frequency_lc))

	start = time.time()

	znow = redshift[isnap]
	znext = redshift[isnap+1]

	#Load the data from Meraxes for this snapshot
	galdata = ReadGalaxiesSnap(isnap,ncores,totngalaxies[isnap],inputdatadir,galdesiredfields)

	# this is the distance needed
	distance = (cosmo.comoving_distance(znow - znext)).value
	
	if galdata["BlackHoleMass"].size>0:

		total_flux = np.zeros(len(galdata["Pos"]))

		position = galdata["Pos"]/cosmo.h
		
		spin_parameter = galdata["Spin"]

		SFR = galdata["Sfr"] ## unit of solar mass/yr

		BHMass = galdata["BlackHoleMass"].astype("float64") #* 1e10 * const.M_sun ## initially 1e10 solar mass but now, in unit of kg
		accretion_rate_hot = galdata["BlackHoleAccretedHotMass"]  / (galdata["dt"] * 1e6 * 365 * 24 * 60 * 60 * un.s)

		accretion_rate_hot /= M_Eddington(BHMass) #galdata["BlackHoleAccretedHotMass"] + galdata["BlackHoleAccretedColdMass"])#
		accretion_rate_hot = accretion_rate_hot.value
		accretion_rate_hot[accretion_rate_hot>1] = 1

		accretion_rate_cold = galdata["BlackHoleAccretedColdMass"] / (galdata["dt"] * 1e6 * 365 * 24 * 60 * 60 * un.s) 
		accretion_rate_cold /= M_Eddington(BHMass) #galdata["BlackHoleAccretedHotMass"] + galdata["BlackHoleAccretedColdMass"])#
		accretion_rate_cold = accretion_rate_cold.value
		accretion_rate_cold[accretion_rate_cold>1] = 1
		
		Lstar_sel = SFR!=0.0
		Lstar = Lnu_SFR(SFR[Lstar_sel])

		flux_star = fluxDensity(Lstar,znow) # this is in Jy at 1.4 GHz

		L_BH_hot_sel = accretion_rate_hot!=0.0
		L_BH_hot = Lnu_BH_radio(BHMass[L_BH_hot_sel] * 1e10 * const.M_sun, const.M_sun, accretion_rate_hot[L_BH_hot_sel], spin_parameter[L_BH_hot_sel])
		

		flux_BH_hot = fluxDensity(L_BH_hot,znow) # this is in Jy at 1.4 GHz

		L_BH_cold_sel = accretion_rate_cold!=0.0
		L_BH_cold = Lnu_BH_quasar(BHMass[L_BH_cold_sel] * 1e10 * const.M_sun, const.M_sun, accretion_rate_cold[L_BH_cold_sel], spin_parameter[L_BH_cold_sel])
		

		flux_BH_cold = fluxDensity(L_BH_cold,znow) # this is in Jy at 1.4 GHz

		total_flux[Lstar_sel] += flux_star.value
		total_flux[L_BH_hot_sel] += flux_BH_hot.value
		total_flux[L_BH_cold_sel] += flux_BH_cold.value

		position = position[total_flux>0]
		total_flux = total_flux[total_flux>0]

		# for now, randomize which direction we take and take the slices corresponding to the distance of lightcone between redshift
		direction = np.random.randint(0,3)
		dimensions = np.arange(0,3,1)[np.arange(0,3,1) != direction]

		# this gives the indices of galaxies that are within the lightcone slice
		ind = indicesRandomization(position[:,direction], bins_pos, distance)

		# select galaxies that are within the slices
		xpos = position[:,dimensions[0]][ind][total_flux[ind]<=max_Jy]
		ypos = position[:,dimensions[1]][ind][total_flux[ind]<=max_Jy]
		fluxpos = total_flux[ind][total_flux[ind]<=max_Jy]

		# scale fg to 150 MHz
		flux_150, gamma = scaledFluxDensity(fluxpos)

		for kk in range(len(frequency_lc)):
			lightcone_fg[:,:, kk] += np.histogram2d(xpos, ypos, bins= np.array([bins_pos,bins_pos]), weights = scaledFluxDensity(flux_150, gamma=gamma, frequency=frequency_lc[kk]))[0]
		

	print("Done snap",isnap,"in",time.time()-start)


# lcfilename = "data/%i/lightcone_galaxies.npz" %random_seed
# np.savez_compressed(lcfilename, lightcone = lightcone)
# lightcone_fg = np.load(lcfilename)["lightcone"]


##########################################################################
#####################   GeneratePowerSpectrum  ###########################
##########################################################################

#Read the lightcone from the meraxes grid data
gridata = ReadGridData(inputdatadir,lcgridfilename,lcdesiredfields)

#Keep only parts that are within the given frequency range
lightcone_eor = gridata["LightconeBox"][:,:,lc_slice_sel]

#The number of frequencies
nfreq = len(frequency_lc)

# the lightcone cube that we fill out with foregrounds
fg_lightcone = np.zeros([Ngrid, Ngrid,nfreq])

# convert foreground to brightness temperature and fill up a cube with its evolution at each frequency
ii = 0
for jj in range(nfreq):
		fg_lightcone[:,:,ii] = lightcone_fg[:,:,jj] / (mpctoRadian(redshift_lc[jj], 500) / Ngrid) / mKtoJy_per_sr(frequency_lc[jj])
		ii+=1

# sort out the coords of box in Mpc
Lz_Mpc = cosmo.comoving_transverse_distance(redshift_lc[-1]).value-cosmo.comoving_transverse_distance(redshift_lc[0]).value + np.diff(cosmo.comoving_transverse_distance(redshift_lc).value)[0]	

coords_xy = np.linspace(0, 500, Ngrid+1)
coords_xy = coords_xy[1:] - np.diff(coords_xy)[0]/2
coords_z = cosmo.comoving_transverse_distance(redshift_lc).value

new_coords_z = np.linspace(coords_z.min(), coords_z.max(), 100)

xx, yy, zz = np.meshgrid(coords_xy, coords_xy, new_coords_z, indexing='ij')

# interpolate the foreground across frequency to smoothen its evolution
f_lc_fg = scipy.interpolate.RegularGridInterpolator([coords_xy, coords_xy, coords_z], fg_lightcone)
highres_lightcone_fg = f_lc(np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T).reshape(Ngrid,Ngrid,len(new_coords_z))

# interpolate the eor across frequency to make sure we have the same dimension as foregrounds
f_lc_fg = scipy.interpolate.RegularGridInterpolator([coords_xy, coords_xy, coords_z], lightcone_eor)
highres_lightcone_eor = f_lc(np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T).reshape(Ngrid,Ngrid,len(new_coords_z))

# Find the PS of foregrounds plus eor
P_both, kperp, kparal = get_power((highres_lightcone_fg + highres_lightcone_eor) * signal.blackmanharris(len(new_coords_z)), [500, 500, Lz_Mpc], bins = 50, res_ndim =2, bin_ave=False, get_variance=False)

# Find the PS of foregrounds only
P_fg = get_power(highres_lightcone_fg * signal.blackmanharris(len(new_coords_z)), [500, 500, Lz_Mpc], bins = 50, res_ndim =2, bin_ave=False, get_variance=False)[0]

# Find the PS of eor only
P_eor = get_power(highres_lightcone_eor * signal.blackmanharris(len(new_coords_z)), [500, 500, Lz_Mpc], bins = 50, res_ndim =2, bin_ave=False, get_variance=False)[0]

#Save the power spectrum
np.savez_compressed(outputdatadir + "power_spectrum.npz", power_spectrum_both = P_both, power_spectrum_eor = P_eor, power_spectrum_fg = P_fg, kperp = kperp, kparal = kparal)
