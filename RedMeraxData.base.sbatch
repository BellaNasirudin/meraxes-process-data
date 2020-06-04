#!/bin/sh
#SBATCH --partition=gstar
#SBATCH -C largemem
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=2:00:00
#SBATCH --account=oz009
#SBATCH --mem=105GB
#SBATCH --output=logs/JOBNAME.log
#SBATCH --error=logs/JOBNAME.error
#SBATCH --job-name=JOBNAME
#SBATCH --tmp=50GB

module load gcc/7.3.0 gsl/2.4 fftw/3.3.7 anaconda3/5.1.0

export PYTHONPATH="/home/bnasirud/PYENV/lib/python3.6/site-packages/"

RUNCMD
