#!/bin/bash

#SBATCH -p normal
#SBATCH -J SW_night2	      	# Job name
#SBATCH -o job.%j.out         	# Name of stdout output file (%j expands to jobId)
#SBATCH -e job.%j.err	      	# Name of error file (%j expands to jobID)
#SBATCH -N 8                  	# Total number of nodes requested
#SBATCH --ntasks-per-node 68 	# Total number of mpi tasks requested (1 node : 68 CPU)
#SBATCH -t 48:00:00          	# Run time (hh:mm:ss), 12  hours
#SBATCH --mail-type=ALL       	# mail alert
#SBATCH --mail-user=yunjaeh@stanford.edu 

echo "Running [your app] with: $cmd on $SLURM_JOB_NODELIST in directory "`pwd`
date

# Working directory
# /scratch/07279/tg865790/BangladeshLES/validationCases/SW_day1
cd $SCRATCH
cd ./BangladeshLES/validationCases/SW_night2

fInput='12_charles.in'
fLog='log.charles12'

ibrun ./solver/charles.exe -i $fInput >  $fLog



