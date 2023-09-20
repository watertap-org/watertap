#!/bin/bash 
#SBATCH --nodes=1  # Run the tasks on the same node
#SBATCH --ntasks-per-node=104 # Tasks per node to be run
#SBATCH --time=1:00:00   # Required, estimate 5 minutes
#SBATCH --account=hpcapps # Required
#SBATCH --partition=debug
#SBATCH --mail-user=kinshuk.panda@nrel.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE

cd /home/kpanda/NAWI/watertap/tutorials/parameter_sweep_demo
module purge
module load craype-x86-spr
module load gcc/13.1.0 anaconda3/2022.05 netlib-lapack/3.11.0-gcc
conda activate /projects/hpcapps/kpanda/conda-envs/watertap

mkdir -p outputs
N_SAMPLES=5000
for nprocs in 100 # {10..100..10}
do
    echo $NSAMPLE $nprocs
    python parameter_sweep_demo_script.py $N_SAMPLES $nprocs > outputs/fout_mp_${N_SAMPLES}_${nprocs} 2> outputs/errout__mp_${N_SAMPLES}_${nprocs}
done
