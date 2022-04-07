#!/bin/bash
#SBATCH -A INF22_npqcd_0
#  #SBATCH -A INF22_test_1
#SBATCH -p m100_usr_prod
#SBATCH --qos=m100_qos_dbg
#SBATCH --error=err_suqa_qms
#SBATCH --output=stdTrDist_log_qms_%j
#SBATCH --time 2:00:00     # format: HH:MM:SS
#SBATCH -N 1                # 1 node
#SBATCH --ntasks-per-node=8 # 8 tasks out of 128
#SBATCH --gres=gpu:4        # 1 gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --exclusive
#SBATCH --job-name=suqa_qms_b0.25_nqe1_e-1.0_3.0


module load gnu/8.4.0
module load blas/3.8.0--gnu--8.4.0
module load lapack/3.9.0--gnu--8.4.0
module load cuda/11.0

#export MPIINC="/cineca/prod/opt/compilers/pgi/19.10/binary/linuxpower/19.10/mpi/openmpi-3.1.3/include/"
#export MPILIB="/cineca/prod/opt/compilers/pgi/19.10/binary/linuxpower/19.10/mpi/openmpi-3.1.3/lib/"

#rm -f stop

#export PGI_ACC_BUFFERSIZE=$((576*16*48*48*(2*2+12)))

beta=0.25
gbeta=1.0
emin=-1.0
emax=3.0
nqe=1
walltime=7200

num_jobs=8
ctr=0
for i in {3..10}
do
  subtag=nlatest_b"$beta"_rt"$i"_nqe"$nqe"_erange"$emin"_"$emax"
  data_fname=rangedef_qms_systest_"$subtag"
	if [ ! -f "$data_fname" ]
	then
		CUDA_VISIBLE_DEVICES=$((ctr%4))  ./qms $beta $gbeta 10000000 "$i" 3 "$nqe" "$data_fname" --max-reverse 1000 --PE-steps 1 --ene-min $emin --ene-max $emax --thermalization 10000 --walltime $walltime --tag "$subtag" >> nlatest_log_qms_stdTrDist_"$subtag" &
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
