num_jobs=4
for i in {1..20}
do
	while (( ${num_jobs@P} >= num_procs )); do
    		wait -n
  	done
    	./qms 0.25 1.0 20000000 "$i" 3 1 qms_systest_"$i" --max-reverse 1000 --PE-steps 1 --ene-min -1 --ene-max 3 --thermalization 10000 &
done
