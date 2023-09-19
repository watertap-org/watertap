#!/bin/sh
mkdir -p outputs
N_SAMPLES=500
for nprocs in {1..10}
do
    # echo $NSAMPLE $nprocs
    python parameter_sweep_demo_script.py $N_SAMPLES $nprocs > outputs/fout_mp_${N_SAMPLES}_${nprocs} 2> outputs/errout__mp_${N_SAMPLES}_${nprocs}
done