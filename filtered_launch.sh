b=0.25; nqe=1; em=-1.0; eM=3.0; for rt in {6..10};do ./qms "$b" 0.8 1000000 $rt 3 "$nqe" test_rdef_nfrutri_b"$b"_gb0.8_em"$em"_eM"$eM"_nqe"$nqe"_rt"$rt" --tag b"$b"_em"$em"_eM"$eM"_nqe"$nqe"_rt"$rt" --max-reverse 100 --ene-min "$em" --ene-max "$eM" --PE-steps 1 --thermalization 10 --walltime 600 &>> log_rt"$rt" && echo "done" && say "Job finished" & done


for nqe in {1..2};do awk FNR!=1 $(ls -v best_systest_test_rdef_nfrutri_b0.25_gb0.8_em-1.1_eM3.1_nqe"$nqe"_rt*) > alldata_b0.25_gb0.8_em-1.1_eM3.1_nqe"$nqe"; done

for nqe in {1..1};do awk FNR!=1 $(ls -v best_systest_test_rdef_nfrutri_b0.25_gb0.8_em-1.0_eM3.0_nqe"$nqe"_rt*) > alldata_b0.25_gb0.8_em-1.0_eM3.0_nqe"$nqe"; done
