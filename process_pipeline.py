from read_meraxes import ReadHeader, ReadGridData
from post_process_parallel import ConstructLightconeFG, GeneratePowerSpectrum 
from conversion import redshifttoFrequency
import numpy as np 
import multiprocessing as mp 
from astropy.cosmology import Planck15 as cosmo
import astropy.units as un
import sys
import os

#The inputdirectory
inputdatadir = sys.argv[1]

#Generate a random number for this run
meraxesdataflag = int(sys.argv[2])

#The output directory
outputdatadir = sys.argv[3]


#Input data
minfrequency = 150e6
maxfrequency = 160e6
startsnap = 60
endsnap = 131
numsnaps = endsnap - startsnap+1
boxsize = 500/cosmo.h*un.Mpc
Ngrid = 256
min_Jy = 1e-8
max_Jy = 10e-3
galdesiredfields = ["dt","Pos","BlackHoleMass","BlackHoleAccretedHotMass","BlackHoleAccretedColdMass","Spin","Sfr"]
lcdesiredfields = ["LightconeBox"]

#Make the output directory with the name in the rseed
if(meraxesdataflag):
	lcgridfilename = "meraxes_grids_48.hdf5"
else:
	lcgridfilename = "meraxes_lightcone.hdf5"

#Read the number of cores and redshift from the main meraxes file
if(meraxesdataflag):
	ncores,redshift,totngalaxies = ReadHeader(inputdatadir,startsnap,endsnap)
else:
	ncores=0
	redshift = np.zeros(numsnaps)
	totngalaxies = np.zeros(numsnaps)

#Read the lightcone from the meraxes grid data
gridata = ReadGridData(inputdatadir,lcgridfilename,["lightcone-z"])

#Get the redshift and frequency of the lightcone
redshift_lc = gridata["lightcone-z"]
frequency_lc = redshifttoFrequency(redshift_lc)
lc_slice_sel=(frequency_lc<=maxfrequency) & (frequency_lc>=minfrequency)
redshift_lc = redshift_lc[lc_slice_sel]
frequency_lc = frequency_lc[lc_slice_sel]

#Generate lightcone from meraxes data
lightcone_fg = ConstructLightconeFG(startsnap,endsnap,inputdatadir, redshift, totngalaxies, ncores,outputdatadir,boxsize,Ngrid,frequency_lc,min_Jy,max_Jy,meraxesdataflag,galdesiredfields)

#Generate the power spectrum
GeneratePowerSpectrum(lightcone_fg,redshift_lc,frequency_lc,lc_slice_sel,Ngrid,inputdatadir,outputdatadir,lcgridfilename,lcdesiredfields)


