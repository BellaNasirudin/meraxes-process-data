import h5py
import numpy as np


def ReadHeader(inputdatadir,startsnap,endsnap):

	#Set the filename
	filename = inputdatadir  + "meraxes.hdf5"

	#Open up file
	allfile = h5py.File(filename,"r")

	#Extract the number of cores meraxes was run with
	ncores = allfile.attrs["NCores"][0]

	# The redshift dataset
	redshift = np.zeros(endsnap+1)

	# The total number of galaxies per snapshot
	totngalaxies = np.zeros(endsnap+1,dtype=np.int32)

	for isnap in range(startsnap,endsnap+1):

		snapkey = "Snap%03d" %isnap

		#Extract the redshift
		redshift[isnap] = allfile[snapkey].attrs["Redshift"][...]

		#Extract the number of galaxies
		totngalaxies[isnap] = allfile[snapkey].attrs["NGalaxies"][...]

	#Close the file
	allfile.close()

	return ncores,redshift,totngalaxies

def ReadGalaxiesSnap(snap,ncores,ngalaxies,inputdatadir,desiredfields):

	#Set the start and end index
	startindx = 0
	endindx = 0

	#Setup the galaxy data
	galdata = {}
	for key in desiredfields:
		if(key=="Pos"):
			galdata[key] = np.zeros([ngalaxies,3])
		else:
			galdata[key] = np.zeros(ngalaxies)

	snapkey = "Snap%03d" %snap

	#Loop over the ncores
	for icore in range(ncores):

		#Set the filename
		filename = inputdatadir + "meraxes_%d.hdf5" %icore

		#Open up the file
		galfile = h5py.File(filename,"r")

		#Update the endindx
		endindx += galfile[snapkey]["Galaxies"].size

		#Load in the desired fields
		for key in desiredfields:
			galdata[key][startindx:endindx] = galfile[snapkey]["Galaxies"][key][:]

		#Update the startind
		startindx = endindx

		#Close the file
		galfile.close()

	return galdata


def ReadGridData(inputdatadir,filename,desiredfields):

	#Open up the file
	lcfile = h5py.File(inputdatadir+filename,"r")

	#Load the data into a dictionary
	lcdata = {}

	#Load in the desiredfielfd
	for key in desiredfields:
		lcdata[key] = lcfile[key][:]

	#Close the file
	lcfile.close()

	return lcdata


def ReadReducedDataSnap(snap,inputdatadir,desiredfields):

	#Open up the file
	snapfile = h5py.File(inputdatadir + "snapshot_%03d.meraxes.hdf5" %snap,"r")

	#Read in the redshift
	redshift = snapfile.attrs["Redshift"][...]
	nextredshift = snapfile.attrs["RedshiftNext"][...]

	#Load the data into a dictionary
	galdata = {}

	#Read in all the data from the file
	for field in desiredfields:
		galdata[field] =  snapfile[field][:]

	return redshift, nextredshift, galdata









