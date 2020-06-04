import numpy as np 
import subprocess
import sys

rand_ints = np.loadtxt("/fred/oz009/bnasirud/meraxes/output/rseed.txt",dtype=int)
script = "reduce_meraxes_data.py"
jobname = "RedMeraxData-"
submitfilename = "RedMeraxData.base.sbatch"
baseoutputdir = "/fred/oz009/bnasirud/reduced_meraxes_data/"
 
#Check if the data has been reduced already
processed_rand_ints = np.loadtxt(baseoutputdir+"rseed.txt" ,dtype=int)


submitfile = open(submitfilename,"r")
submitfiletxt = submitfile.read()
submitfile.close()

for random_seed in rand_ints:

	#Check to see if this random int has already been processed
	if(random_seed in processed_rand_ints):
		print(random_seed,"has already been reduced, skipping...")
		continue

	#Replace the RANDN in the submit file
	newsubmitfiletxt = submitfiletxt.replace("JOBNAME",jobname+str(random_seed))
	newsubmitfiletxt = newsubmitfiletxt.replace("RUNCMD","python3 " + script + " " + str(random_seed) + " " + baseoutputdir)

	#Open up the new submit files
	newsubmitfilename = "runscripts/" + jobname +str(random_seed)+".sbatch"
	newsubmitfile = open(newsubmitfilename,"w")
	newsubmitfile.write(newsubmitfiletxt)
	newsubmitfile.close()

	#Now submit Meraxes using subprocess
	returncode = subprocess.call("sbatch " + newsubmitfilename, shell=True)

	if(returncode!=0):
		raise SystemExit("Submission failed")
	else:
		print(script,"has been submited for random seed",random_seed)
