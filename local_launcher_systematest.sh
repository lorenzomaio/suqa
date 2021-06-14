#num_jobs=1
for i in {1..15}
do
#	while (( ${num_jobs@P} >= num_procs )); do
#    		wait -n
#  	done
        echo $i
    	./qms 0.25 1.0 20000000 "$i" 3 1 qms_systest_b0.25_rt"$i" --max-reverse 1000 --PE-steps 1 --ene-min -1 --ene-max 3 --thermalization 10000 >> log_qms_systest_b0.25_rt$i
done
