from reduce_meraxes_data_func import ReduceMeraxesData
import sys
import os


#The random seed to run with
random_seed = int(sys.argv[1])

#The base output directory
baseoutputdir = sys.argv[2]


#Set the output directory name
outputdir=baseoutputdir+str(random_seed)+"-EFN0.0001/"

startsnap=60
endsnap=131

#The meraxes base directory
meraxesbasedir = "/fred/oz009/bnasirud/meraxes/"

#Set the input and output directories
meraxesoutputdir = "/fred/oz009/bnasirud/meraxes/output/output-Genesis_L500_N2160-"+str(random_seed)+"-EFN0.0/"

#The desired fields from the meraxes files
galdesiredfields = ["Pos","Spin","Sfr","BlackHoleMass","BlackHoleAccretedHotMass","BlackHoleAccretedColdMass","dt"]
griddesiredfields = ["LightconeBox","lightcone-z"]

#Make the directory to put the reduced data
os.makedirs(outputdir,exist_ok=True)

#Run the function to reduced the meraxes data for the desired snapshots
ReduceMeraxesData(startsnap,endsnap,meraxesoutputdir,outputdir,galdesiredfields,griddesiredfields)

#Write seed to file containing the list of random numbers
if(os.path.exists(baseoutputdir + "rseed.txt")):
	mode = "a"
else:
	mode = "w"

randnfile = open(baseoutputdir + "rseed.txt",mode)
randnfile.write(str(random_seed)+"\n")
randnfile.close()