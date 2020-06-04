import h5py
import numpy as np
import shutil

rand_ints = np.loadtxt("/fred/oz009/bnasirud/meraxes/output/rseed.txt",dtype=int)

base_meraxes_file = "/fred/oz009/bnasirud/meraxes/output/output-Genesis_L500_N2160-1809/meraxes.hdf5" 

#The meraxes base directory
meraxesbasedir = "/fred/oz009/bnasirud/meraxes/"

startsnap = 60
endsnap = 132
irmdir = False
deleteindexes = []

for i,random_seed in enumerate(rand_ints):

	#Make the output directory with the name in the rseed
	inputdatadir =meraxesbasedir + "output/" + "output-Genesis_L500_N2160-" + str(random_seed)+"/"

	try:
		hdffile = h5py.File(inputdatadir + "meraxes_0.hdf5","r")
		
		for isnap in range(startsnap,endsnap):
			snapkey = "Snap%03d" %isnap
			try:
				galdata = hdffile[snapkey]["Galaxies"]
			except KeyError:
				print(snapkey,"Does not exist in ",inputdatadir + "meraxes_0.hdf5")
				print("Deleting the folder",inputdatadir)
				shutil.rmtree(inputdatadir)
				irmdir=True
				deleteindexes.append(i)
				break

		hdffile.close()
	except IOError:
		print("Cannot open meraxes_0.hdf5 for seed",random_seed)

	if(irmdir==False):
		try:
			hdffile = h5py.File(inputdatadir + "meraxes.hdf5","r")
			hdffile.close()
		except IOError:
			print("Cannot open meraxes.hdf5 for seed",random_seed)
			shutil.copyfile(base_meraxes_file,inputdatadir + "meraxes.hdf5")
			print("Done copying the base meraxes file")

print("Removing the random_seed from the file",rand_ints[deleteindexes])
rand_ints = np.delete(rand_ints,deleteindexes)

np.savetxt("/fred/oz009/bnasirud/meraxes/output/rseed.txt",rand_ints,fmt="%d")



