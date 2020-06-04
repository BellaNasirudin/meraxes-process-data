import numpy as np 
import subprocess
import sys
import argparse


# The number of meraxes to submit
nsubmit = int(sys.argv[1])

# The escape fraction to run with
escape_fraction_norm = sys.argv[2]


script = "meraxes_pipeline.py"
jobname = "MeraxesPipe-"
submitfilename = "MeraxesPipe.base.sbatch"
meraxesbasedir = "/fred/oz009/bnasirud/meraxes/"
basemeraxesoutputdir = meraxesbasedir + "output/"
basereducedoutputdir = "/fred/oz009/bnasirud/reduced_meraxes_data/"
baseprocessedoutputdir = "/fred/oz009/bnasirud/data_processed/"
submitcmd = "sbatch"

submitfile = open(submitfilename,"r")
submitfiletxt = submitfile.read()
submitfile.close()

for i in range(nsubmit):

	#Get a random number to set the Meraxes seed
	random_seed = np.random.randint(0,int(2**32-1))

	meraxesoutputdir = basemeraxesoutputdir + "output-Genesis_L500_N2160-"+ str(random_seed) +"-EFN" + escape_fraction_norm + "/"
	reducedoutputdir = basereducedoutputdir + str(random_seed) +"-EFN" + escape_fraction_norm + "/"
	processedoutputdir  = baseprocessedoutputdir + str(random_seed) +"-EFN" + escape_fraction_norm + "/"

	#Replace the RANDN in the submit file
	newsubmitfiletxt = submitfiletxt.replace("JOBNAME",jobname+str(random_seed))
	newsubmitfiletxt = newsubmitfiletxt.replace("RUNCMD","python3 " + script +" " + meraxesbasedir + " " + str(random_seed) + " " + escape_fraction_norm + " " + meraxesoutputdir + " " + reducedoutputdir + " " + processedoutputdir)

	#Open up the new submit files
	newsubmitfilename = "runscripts/"+jobname+str(random_seed)+".sbatch"
	newsubmitfile = open(newsubmitfilename,"w")
	newsubmitfile.write(newsubmitfiletxt)
	newsubmitfile.close()

	#Now submit Meraxes using subprocess
	returncode = subprocess.call(submitcmd + " " + newsubmitfilename, shell=True)

	if(returncode!=0):
		raise SystemExit("Submission failed")
	else:
		print(script,"has been submited for random seed",random_seed,"and EFN",escape_fraction_norm)

	np.random.seed(random_seed)