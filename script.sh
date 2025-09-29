#!/bin/bash
#SBATCH --job-name=PR10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
##SBATCH --mem=0
#SBATCH --mem-per-cpu=10G
##SBATCH --time=0-12:00:00
#SBATCH --partition=serial#highmem
##SBATCH --array=0-3

module load julia/1.5.0

#export Rstart=3.4
#export dR=0.2
#export N=1
#
#export DIST=$(echo $Rstart-$SLURM_ARRAY_TASK_ID*$dR |bc)

#export Ort=$(pwd)
#export EXEC=/home/tserwatk/DMRG/MBpol/dE/
#
##export NEWDIR=$Ort/tmp_$DIST
##mkdir -p $NEWDIR
#cp $EXEC/*.jl $Ort
#cp $EXEC/test-mbpol $Ort

#sed "s/startValue/$DIST/;s/stepSize/$dR/;s/numberOfCalculations/$N/" $Ort/input.txt > $NEWDIR/input.txt

#cd $NEWDIR

julia main.jl

