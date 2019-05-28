#!/bin/bash

#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=50:00:00
#PBS -q batch
#PBS -j oe
#PBS -o /storage/users/fbaharifard/ComplexNetworks/out.out
#PBS -N Runme 
#PBS -m ea
#PBS -M fateme.baharifard@gmail.com 


export PATH=/share/Application/anaconda3/bin:$PATH
export LD_LIBRARY_PATH=/share/Application/anaconda3/lib:$LD_LIBRARY_PATH

cd $PBS_O_WORKDIR

OMP_NUM_THREADS=$PBS_NUM_PPN
export  OMP_NUM_THREADS
echo 'numberOfThreads=' $OMP_NUM_THREADS

#-------------------User Define ------------------
export InputName=runme

#python runme.py
python ContinueRunme.py



