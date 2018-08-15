#!/bin/bash

#SBATCH --workdir=./ # Working directory
#SBATCH --mail-user=battistn@tcnj.edu       # Who to send emails to
#SBATCH --mail-type=ALL             # Send emails on start, end and failure
#SBATCH --job-name=matlab	    # Job name
#SBATCH --output=job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH --nodes=1                   # Total number of nodes (a.k.a. servers) requested
#SBATCH --ntasks=1                  # Total number of mpi tasks requested
#SBATCH --partition=nolimit         # Partition (a.k.a.queue) to use
#SBATCH --time=10-00:00:00          # Run time (days-hh:mm:ss)
#SBATCH --nodelist=node006 

# Launch a serial job
echo "Starting @ "`date`
matlab -nodisplay < main2d.m > main2d.out
echo "Completed @ "`date`
