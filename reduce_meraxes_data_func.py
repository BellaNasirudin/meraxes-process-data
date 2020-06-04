import h5py
import time
import sys
import os
import shutil
from read_meraxes import ReadHeader, ReadGalaxiesSnap
import multiprocessing as mp
from itertools import repeat
import numpy as np 

# os.system("taskset -p 0xff %d" % os.getpid())


def WriteSnapshot(isnap,endsnap,ngalaxies,redshift,inputdatadir,ncores,outfilebasename,desiredfields):

	start = time.time()

	#Open up the file to output the data to
	outfilename = outfilebasename %isnap
	outfile = h5py.File(outfilename,"w")

	#Write the redshift and number of galaxies for the redshift
	
	outfile.attrs["Redshift"] = redshift[isnap]
	if(isnap!=endsnap): outfile.attrs["RedshiftNext"] = redshift[isnap+1]
	outfile.attrs["NGalaxies"] = ngalaxies


	galdata = ReadGalaxiesSnap(isnap,ncores,ngalaxies,inputdatadir,desiredfields)

	sel = (galdata["Sfr"][:]>0) | (galdata["BlackHoleMass"][:]>0)

	# Write all the desired fields
	for field in desiredfields:
		outfile.create_dataset(field,data=galdata[field][sel],compression="gzip",compression_opts=6)

	outfile.close()

	print("Done snap",isnap,"in",time.time()-start,flush=True)

def ReduceMeraxesData(startsnap,endsnap,inputdatadir,outputdir,galdesiredfields,griddesiredfields):

	#Input file names
	allfilename = inputdatadir + "meraxes.hdf5"
	gridfilename = inputdatadir + "meraxes_grids_48.hdf5"

	#Outputfilenames
	outfilebasename = outputdir + "snapshot_%03d.meraxes.hdf5"
	lightconefilename = outputdir + "meraxes_lightcone.hdf5"

	totstart = time.time()

	#Read the number of cores and redshift from the main meraxes file
	ncores,redshift,totngalaxies = ReadHeader(inputdatadir,startsnap,endsnap)

	#The snapshot to run over
	snapshots = np.arange(startsnap,endsnap+1)

	#Get the number of cpus in this job if the varible exists otherwise set it as the total number of cpus
	ncpus = int(os.getenv("SLURM_NTASKS",mp.cpu_count()))

	#Setup a pool to create
	pool = mp.Pool(ncpus)

	#Run the pool to process the snapshot
	pool.starmap(WriteSnapshot,zip(snapshots,repeat(endsnap),totngalaxies[snapshots],repeat(redshift),repeat(inputdatadir),repeat(ncores),repeat(outfilebasename),repeat(galdesiredfields)))

	### Extract the data from the lightconefile ###

	#Read the lightcone and its redshift
	gridfile = h5py.File(gridfilename,"r")

	#Open up a file to store the lightcone data
	lightconefile = h5py.File(lightconefilename,"w")

	#Extract the data and write it to a new file
	for key in griddesiredfields:
		lightconefile.create_dataset(key,data=gridfile[key][:],compression="gzip",compression_opts=6)

	gridfile.close()
	lightconefile.close()

	shutil.copyfile(inputdatadir+"/Genesis-L500_N2160-parameters.par",outputdir+"/Genesis-L500_N2160-parameters.par")
	shutil.copyfile(inputdatadir+"/Genesis-L500_N2160-input.par",outputdir+"/Genesis-L500_N2160-input.par")

	print("Done reducing the data in",time.time()-totstart)
