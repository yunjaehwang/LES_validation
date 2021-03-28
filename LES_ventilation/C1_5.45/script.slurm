#!/bin/bash

#SBATCH -p normal
#SBATCH -J C1_5.45	      	# Job name
#SBATCH -o job.%j.out         	# Name of stdout output file (%j expands to jobId)
#SBATCH -e job.%j.err	      	# Name of error file (%j expands to jobID)
#SBATCH -N 2                  	# Total number of nodes requested
#SBATCH --ntasks-per-node 68 	# Total number of mpi tasks requested (1 node : 68 CPU)
#SBATCH -t 48:00:00          	# Run time (hh:mm:ss), 12  hours
#SBATCH --mail-type=ALL       	# mail alert
#SBATCH --mail-user=yunjaeh@stanford.edu 

echo "Running [your app] with: $cmd on $SLURM_JOB_NODELIST in directory "`pwd`
date

# Working directory
# /scratch/07279/tg865790/LES_validation/LES_ventilation/E1_10
cd $SCRATCH
cd ./LES_validation/LES_ventilation/C1_5.45

fInput='04_charles.in'
fLog='log.charles04'

ibrun $HOME/CharLES/nextgen_turb/charles.exe -i $fInput > $fLog

