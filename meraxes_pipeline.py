from execute_meraxes import run_meraxes
from reduce_meraxes_data_func import ReduceMeraxesData
from read_meraxes import ReadHeader, ReadGridData
from post_process_parallel import ConstructLightconeFG, GeneratePowerSpectrum 
from conversion import redshifttoFrequency
import sys
import os
from astropy.cosmology import Planck15 as cosmo
import astropy.units as un
import numpy as np

#The meraxes base directory
meraxesbasedir = sys.argv[1]

#The random seed to run with
random_seed = sys.argv[2]

#The escape fraction to run with
escape_fraction_norm = sys.argv[3]

#The meraxes output directory 
meraxesoutputdir = sys.argv[4]

#The reduced data output directory
reducedoutputdir = sys.argv[5]

#The processed data output directory
processedoutputdir = sys.argv[6]

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
lcgridfilename = "meraxes_lightcone.hdf5"
galdesiredfields = ["dt","Pos","BlackHoleMass","BlackHoleAccretedHotMass","BlackHoleAccretedColdMass","Spin","Sfr"]
lcdesiredfields = ["LightconeBox","lightcone-z"]

# ################################################
# ############### Running meraxes ################
# ################################################


#Make the directory for the meraxes data to be put
os.makedirs(meraxesoutputdir,exist_ok=True)

#First run meraxes
returncode = run_meraxes(meraxesbasedir,meraxesoutputdir,random_seed,escape_fraction_norm)

#Check if meraxes completed successfully
if(returncode!=0):
	raise SystemExit("Running meraxes failed")

#Write seed to file containing the list of random numbers
rseedname=meraxesbasedir + "output/rseed-EFN" + escape_fraction_norm +".txt"
if(os.path.exists(rseedname)):
	mode = "a"
else:
	mode = "w"

randnfile = open(rseedname,mode)
randnfile.write(random_seed+"\n")
randnfile.close()


#################################################
########### Reducing the meraxes data ###########
#################################################

#Make the directory to put the reduced data
os.makedirs(reducedoutputdir,exist_ok=True)

#Run the function to reduced the meraxes data for the desired snapshots
ReduceMeraxesData(startsnap,endsnap,meraxesoutputdir,reducedoutputdir,galdesiredfields,lcdesiredfields)

#Write seed to file containing the list of random numbers
rseedname ="/fred/oz009/bnasirud/reduced_meraxes_data/rseed-EFN" + escape_fraction_norm +".txt"
if(os.path.exists(rseedname)):
	mode = "a"
else:
	mode = "w"

randnfile = open(rseedname,mode)
randnfile.write(random_seed+"\n")
randnfile.close()

#################################################
############## Processing Pipeline ##############
#################################################

#Make the directory to put the processed data
os.makedirs(processedoutputdir,exist_ok=True)

#Read the data from the Reduced data
meraxesdataflag=0
ncores=0
redshift = np.zeros(numsnaps)
totngalaxies = np.zeros(numsnaps)

#Read the lightcone from the meraxes grid data
gridata = ReadGridData(reducedoutputdir,lcgridfilename,["lightcone-z"])

#Get the redshift and frequency of the lightcone
redshift_lc = gridata["lightcone-z"]
frequency_lc = redshifttoFrequency(redshift_lc)
lc_slice_sel=(frequency_lc<=maxfrequency) & (frequency_lc>=minfrequency)
redshift_lc = redshift_lc[lc_slice_sel]
frequency_lc = frequency_lc[lc_slice_sel]

#Make the directory to put the processed data
os.makedirs(processedoutputdir,exist_ok=True)

#Generate lightcone from meraxes data
lightcone_fg = ConstructLightconeFG(startsnap,endsnap,reducedoutputdir, redshift, totngalaxies, ncores,processedoutputdir,boxsize,Ngrid,frequency_lc,min_Jy,max_Jy,meraxesdataflag,galdesiredfields)

#Generate the power spectrum
GeneratePowerSpectrum(lightcone_fg,redshift_lc,frequency_lc,lc_slice_sel,Ngrid,reducedoutputdir,processedoutputdir,lcgridfilename,lcdesiredfields)