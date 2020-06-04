import h5py
import numpy as np
import astropy.constants as const
from astropy.cosmology import Planck15 as cosmo
import astropy.units as un

# Radio luminosity of BH
# Based on Amarantidis et al (2019): https://academic.oup.com/mnras/article/485/2/2694/5365431
def Lnu_BH_radio(BHMass, SolarMass, accretion_rate_hot, spin_parameter, nu=1.4e9): # needs hot gas!

	if(BHMass.size==1):

		normalization = 8e-5

		Ljet_radiomode = 2 * 1e45 * BHMass / 1e9 / SolarMass * accretion_rate_hot / 0.01 * spin_parameter**2

		Lnu = normalization * (BHMass / 1e9 / SolarMass * accretion_rate_hot / 0.01)**(0.42) * Ljet_radiomode / nu

		return Lnu.value # W/Hz
	else:

		Lnu = np.zeros(BHMass.size)

		#### Radio Mode ####

		sel =  accretion_rate_hot !=0

		normalization = 8e-5


		Ljet_radiomode = 2 * 1e45 * BHMass[sel] / 1e9 / SolarMass * accretion_rate_hot[sel] / 0.01 * spin_parameter[sel]**2


		Lnu[sel] = normalization * (BHMass[sel] / 1e9 / SolarMass * accretion_rate_hot[sel] / 0.01)**(0.42) * Ljet_radiomode / nu


		return Lnu

def Lnu_BH_quasar(BHMass, SolarMass, accretion_rate_cold, spin_parameter, nu=1.4e9):

	if(BHMass.size==1):

		normalization = 5e-2

		Ljet_quasar = 2.5 * 1e43 * (BHMass / 1e9 / SolarMass)**(1.1) * (accretion_rate_cold / 0.01)**(1.2) * spin_parameter**2

		Lnu = normalization * (BHMass / 1e9 / SolarMass)**(0.32) * (accretion_rate_cold / 0.01)**(-1.2) * Ljet_quasar / nu

		return Lnu.value # W/Hz

	else:

		Lnu = np.zeros(BHMass.size)

		sel =  accretion_rate_cold !=0

		normalization = 5e-2

		Ljet_quasar = 2.5 * 1e43 * (BHMass[sel] / 1e9 / SolarMass)**(1.1) * (accretion_rate_cold[sel] / 0.01)**(1.2) * spin_parameter[sel]**2


		Lnu[sel] = normalization * (BHMass[sel] / 1e9 / SolarMass)**(0.32) * (accretion_rate_cold[sel] / 0.01)**(-1.2) * Ljet_quasar / nu

		return Lnu


# Conversion of mass to Eddington Luminosity
# Based on Fanidakis et. al (2011) https://academic.oup.com/mnras/article/410/1/53/1031503
def M_Eddington(BHMass, accretionEfficiency=0.1):

	L_Eddington = 4 * np.pi * const.G * BHMass * const.c  * const.m_p / const.sigma_T

	rate_Eddington = L_Eddington / const.c**2 / accretionEfficiency
	
	return rate_Eddington # 1e10 solar mass /s

# Flux density in Jansky
def fluxDensity(Luminosity, redshift):

	distance = cosmo.luminosity_distance(redshift).to('m')

	prefactor = 1 / (4 * np.pi *  distance**2 ) / (1+redshift)**(0.8-1)
	
	return 1e26 * prefactor * Luminosity # 10^26 W/m^2/Hz = 1 Jy

# From Condon 1992 
# Gives luminosity at 1.4 GHz
def Lnu_SFR(SFR, nu_GHz=1.4):

	return SFR * (5.5e20 * (nu_GHz)**(-0.1) + 5.3e21 * (nu_GHz)**(-0.8)) #from condon 1992
	#/ (7.5 * 10**(-29)) # erg /s /Hz https://www.aanda.org/articles/aa/pdf/2002/35/aa2265.pdf
	# * 8.07 * 1e20 # unit of W/Hz based on https://academic.oup.com/mnras/article/308/1/45/1007213

def redshifttoFrequency(redshift):
	Frequency = 1420e6 / (redshift + 1)
	
	return Frequency

def frequencytoRedshift(frequency):
	redshift = 1420e6 / frequency - 1
	
	return redshift

def scaledFluxDensity(fluxDensity, frequency=1.4e9, gamma= None, redshift = None, frequency0 = 150e6, mu=0.8, sigma=0.05):

	if redshift != None:
		frequency0 = frequencyRedshift(redshift)

	if gamma is None:
		gamma = np.random.normal(mu,sigma, len(fluxDensity))
		
		return fluxDensity * (frequency/frequency0)**(-gamma), gamma #in Jy

	else:

		return fluxDensity * (frequency/frequency0)**(-gamma) #in Jy

def mpctoRadian(boxsize, redshift):
	return boxsize / cosmo.comoving_transverse_distance(redshift).value


def mKtoJy_per_sr(nu):
	wvlngth = const.c / (nu / un.s)

	intensity = 2 * const.k_B * 1e-3 * un.K / wvlngth ** 2

	flux_density = 1e26 * intensity.to(un.W / (un.Hz * un.m ** 2))

	return flux_density.value



