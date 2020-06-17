import numpy as np 
import subprocess
import sys
import time
import os

escape_fraction_norm = sys.argv[1]

rseedfname = "/fred/oz009/bnasirud/reduced_meraxes_data/rseed-EFN"+escape_fraction_norm+".txt"
# rseedfname = "/fred/oz009/bnasirud/processing_pipeline/reproccess-rseed-EFN0.06.txt"

rand_ints = np.loadtxt(rseedfname,dtype=int)
script = "process_pipeline.py"
jobname = "ProcessPipeline-"
submitfilename = "ProcessPipeline.base.sbatch"
meraxesdataflag = 0 # 0 for reduced data, 1 for meraxes data 

#The meraxes base directory
meraxesbasedir = "/fred/oz009/bnasirud/meraxes/"

submitfile = open(submitfilename,"r")
submitfiletxt = submitfile.read()
submitfile.close()

for random_seed in rand_ints:

	# random_seed = 4125055411 #986284303 #1176542000 #
	#The output directory
	outputdatadir = "/fred/oz009/bnasirud/data_processed/"+str(random_seed)+"-EFN" + escape_fraction_norm + "/" 
	# if(os.path.isdir(outputdatadir)):

	# if(os.path.exists(outputdatadir + "power_spectrum.npz")):
	# 	checkfile = np.load(outputdatadir + "power_spectrum.npz")
	# 	if("power_spectrum_both" in list(checkfile)):
	# 		print("power_spectrum_both already exists in",outputdatadir+ "power_spectrum.npz, skipping...")
	# 		continue
	# 	else:
	# 		print("power_spectrum_both does not exist in",outputdatadir+ "power_spectrum.npz, reprocessing")
	# else:
	# 	os.makedirs(outputdatadir,exist_ok=True)

	#Make the output directory with the name in the rseed
	if(meraxesdataflag):
		inputdatadir = meraxesbasedir + "output/" + "output-Genesis_L500_N2160-" + str(random_seed)+"-EFN" + escape_fraction_norm + "/"
	else:
		inputdatadir = "/fred/oz009/bnasirud/reduced_meraxes_data/" + str(random_seed)+"-EFN" + escape_fraction_norm + "/"

	#Replace the RANDN in the submit file
	newsubmitfiletxt = submitfiletxt.replace("JOBNAME",jobname+str(random_seed))
	newsubmitfiletxt = newsubmitfiletxt.replace("RUNCMD","python3 " + script + " " + inputdatadir + " " + str(meraxesdataflag) + " " + outputdatadir)

	#Open up the new submit files
	newsubmitfilename = "runscripts/" + jobname + str(random_seed) + ".sbatch"
	newsubmitfile = open(newsubmitfilename,"w")
	newsubmitfile.write(newsubmitfiletxt)
	newsubmitfile.close()

	#Start the timer
	start = time.time()

	#Now submit Meraxes using subprocess
	returncode = subprocess.call("sbatch " + newsubmitfilename, shell=True)

	if(returncode!=0):
		raise SystemExit("Submission failed")
	else:
		print("Submitted script",script,"for random seed",random_seed,"and EFN",escape_fraction_norm)

	break
