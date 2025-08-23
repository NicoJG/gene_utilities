#!/bin/bash -l
#SBATCH -J GENE              ### job name
#SBATCH --nodes=2            ### Total number of nodes
#SBATCH --ntasks-per-node=64 ### MPI tasks per node (64/node for bonobo, 16/node for edge)
##SBATCH --ntasks=64         ### MPI tasks (only needed on edge)
#SBATCH --cpus-per-task=1    ### Number of threads per task (OpenMP)
#SBATCH --time=24:00:00      ### wall clock time
#SBATCH --mem-per-cpu=16000M  ### memory (bonobo: 16000M, edge: 4000M, viper: 4G)
##SBATCH --partition=dbg ### Uncomment for small jobs
#SBATCH -o ./output/job_%x.%j.out
#SBATCH -e ./output/job_%x.%j.err

# node/cpu setups on different clusters:
# bonobo: 64cpus/node, min 2 nodes -> min 128 tasks
# edge:   16cpus/node
# viper:  128cpus/node

export MEMORY_PER_CORE=16000 # bonobo: 16000, edge: 4000, viper: 4000

## Run this as if it would be in a gene problem directory
PARAMETERS_FILE="./parameters"

# Create temporary directory with job ID to copy the problem template (with symlinks)
TEMPORARY_DIR="$GENEDIR/prob_job_$SLURM_JOB_ID"
TEMPLATE_DIR="$GENEDIR/prob01"

# Function to cleanup temporary directory
cleanup() {
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))
    SECONDS=$((DURATION % 60))
    
    echo "Job ended at: $(date +'%Y-%m-%d %H:%M:%S (UTC%:z)')"
    echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
    echo "Cleaning up temporary directory: $TEMPORARY_DIR"
    if [ -d "$TEMPORARY_DIR" ]; then
        rm -rf "$TEMPORARY_DIR"
    fi
}

# Set trap to cleanup on exit (both success and failure)
trap cleanup EXIT

# Create temporary directory by copying prob_template
cp -r "$TEMPLATE_DIR" "$TEMPORARY_DIR"

# Copy parameters file to temporary directory
cp "$PARAMETERS_FILE" "$TEMPORARY_DIR/parameters"

# Change to temporary directory
cd "$TEMPORARY_DIR"

# Output start time and capture timestamp
START_TIME=$(date +%s)
echo "Job started at: $(date +'%Y-%m-%d %H:%M:%S (UTC%:z)')"

# Default script:

module purge
MODULELIST=$(make -f ../makefile USE_GPU=no get_required_modules | tail -n 1)
eval $MODULELIST
#module list

#FUTILS requirement:
export HDF5_USE_FILE_LOCKING=FALSE

## execute one of the following

MACHINENAME=$(make -f ../makefile USE_GPU=no mach_wrapper)

## single run:
mpirun ./gene_$MACHINENAME

## scan (see GENE documentation or ./scanscript --help):
#./scanscript --np $SLURM_NTASKS --ppn 128 --mps 4 --syscall="mpirun ./gene_$MACHINENAME"

