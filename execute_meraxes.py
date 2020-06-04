import numpy as np 
import subprocess
import sys
import os


def run_meraxes(meraxesbasedir,outputdir,rseed,escape_fraction_norm):

	#The directories
	meraxesexec = meraxesbasedir + "src/build/bin/meraxes"
	meraxesinputdir=meraxesbasedir+"input/"
	meraxestestdir=meraxesbasedir+"test_run/"

	#Input file name
	parameterfilename = "Genesis-L500_N2160-parameters.base.par"
	inputfilename = "Genesis-L500_N2160-input.base.par"

	#Open up the parameter and input files and read their contents
	parameterfile = open(meraxesinputdir+parameterfilename,"r")
	parameterfiletxt = parameterfile.read()
	parameterfile.close()
	inputfile = open(meraxestestdir+inputfilename,"r")
	inputfiletxt = inputfile.read()
	inputfile.close()
	

	#Set the Meraxes rseed in the parameterfile
	newparameterfiletxt = parameterfiletxt.replace("RSEED",rseed)
	newparameterfiletxt = parameterfiletxt.replace("ESCAPEFRACTIONNORM",escape_fraction_norm)

	#Write this to a new paramterfile in the output directory
	newparameterfilename = outputdir+"/Genesis-L500_N2160-parameters.par"
	newparameterfile = open(newparameterfilename,"w")
	newparameterfile.write(newparameterfiletxt)
	newparameterfile.close()

	#Update the directories in the input file
	newinputfiletxt = inputfiletxt.replace("MERAXESBASEDIR",meraxesbasedir)
	newinputfiletxt = newinputfiletxt.replace("OUTPUTDIR",outputdir)
	newinputfiletxt = newinputfiletxt.replace("DEFAULTSFILE",newparameterfilename)

	#Write this to a new inputfile in the output directory
	newinputfilename = outputdir+"/Genesis-L500_N2160-input.par"
	newinputfile = open(newinputfilename,"w")
	newinputfile.write(newinputfiletxt)
	newinputfile.close()

	#Now run Meraxes using subprocess
	return subprocess.call("srun %s %s" %(meraxesexec,newinputfilename)  , shell=True)

	#Output file name
	# randnfilename = "rseed.txt"

	#Write seed to file containing the list of random numbers
	# randnfile = open(meraxesoutputdir + randnfilename,"a")
	# randnfile.write(str(rseed)+"\n")
	# randnfile.close()

















