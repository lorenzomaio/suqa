#!/bin/bash
#SBATCH -A INF21_npqcd_0
#  #SBATCH -A INF20_test_1
#SBATCH -p m100_usr_prod
#SBATCH --qos=m100_qos_dbg
#SBATCH --error=err_suqa_qms
#SBATCH --output=stdTrDist_log_qms_%j
#SBATCH --time 02:00:00     # format: HH:MM:SS
#SBATCH -N 1                # 1 node
#SBATCH --ntasks-per-node=8 # 8 tasks out of 128
#SBATCH --gres=gpu:4        # 1 gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --exclusive
#SBATCH --job-name=suqa_qms_b0.25


module load gnu/8.4.0
module load blas/3.8.0--gnu--8.4.0
module load lapack/3.9.0--gnu--8.4.0
module load cuda/10.1

#export MPIINC="/cineca/prod/opt/compilers/pgi/19.10/binary/linuxpower/19.10/mpi/openmpi-3.1.3/include/"
#export MPILIB="/cineca/prod/opt/compilers/pgi/19.10/binary/linuxpower/19.10/mpi/openmpi-3.1.3/lib/"

#rm -f stop

#export PGI_ACC_BUFFERSIZE=$((576*16*48*48*(2*2+12)))

num_jobs=8
ctr=0
for i in {1..20}
do
	if [ ! -f "qms_systest_$i" ]
	then
		CUDA_VISIBLE_DEVICES=$((ctr%4))  ./qms 0.25 1.0 20000000 "$i" 3 1 qms_systest_"$i" --max-reverse 1000 --PE-steps 1 --ene-min -1 --ene-max 3 --thermalization 10000 --walltime 7100 >> log_qms_stdTrDist_"$i" &
		((ctr++))
	fi

	while :
	do
		background=( $(jobs -p) )
		if (( ${#background[@]} < num_jobs ))
		then 
			break
		fi
		sleep 1
	done
done
wait
