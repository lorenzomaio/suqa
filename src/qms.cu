#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <cstring>
#include <stdio.h>
//#include <bits/stdc++.h>
//#include <unistd.h>
#include <cmath>
#include <cassert>
#include <lapack.h>
#include "Rand.hpp"
#include <chrono>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_device_runtime_api.h>
#include "io.hpp"
#include "parser.hpp"
#include "suqa.cuh"
#include "system.cuh"
#include "qms.cuh"

//XXX: test only, remove after
extern double *host_state_re, *host_state_im;
const double eev[8] = {-1.,-1.,-1.,-1.,-1.,-1.,3.,3.};


const int Amat_coo[16][2] = {{0,6},{0, 7},{1,6},{1,7},{2,4},{2, 5},{3,4},{3,5},{4,2},{4, 3},{5,2},{5,3},{6,0},{6, 1},{7,0},{7,1}};
const double Amat_val[16][2] = {{1,0},{0,-1},{0,1},{1,0},{1,0},{0,-1},{0,1},{1,0},{1,0},{0,-1},{0,1},{1,0},{1,0},{0,-1},{0,1},{1,0}};

using namespace std;





#define NUM_THREADS 128
#define MAXBLOCKS 65535
uint suqa::threads;
uint suqa::blocks;
cudaStream_t suqa::stream1, suqa::stream2;


// simulation parameters
double beta;
double h;
int thermalization;

// defined in src/system.cu
void init_state(ComplexVec& state, uint Dim);

arg_list args;

void save_measures(string outfilename){
    FILE * fil = fopen(outfilename.c_str(), "a");
    for(uint ei = 0; ei < qms::E_measures.size(); ++ei){
        fprintf(fil, "%.16lg %.16lg\n", qms::E_measures[ei], qms::X_measures[ei]);
    }
    fclose(fil);
    qms::E_measures.clear();
    qms::X_measures.clear();
}

void deallocate_state(ComplexVec& state){
    if(state.data!=nullptr){
        HANDLE_CUDACALL(cudaFree(state.data));
    }
    state.vecsize=0U;
}

void allocate_state(ComplexVec& state, uint Dim){
    if(state.data!=nullptr or Dim!=state.vecsize)
        deallocate_state(state);


    state.vecsize = Dim; 
    HANDLE_CUDACALL(cudaMalloc((void**)&(state.data), 2*state.vecsize*sizeof(double)));
    // allocate both using re as offset, and im as access pointer.
    state.data_re = state.data;
    state.data_im = state.data_re + state.vecsize;
}


int main(int argc, char** argv){
    if(argc < 8){
        printf("usage: %s <beta> <g_beta> <metro steps> <reset each> <num state qbits> <num ene qbits> <output file path> [--max-reverse <max reverse attempts> (20)] [--seed <seed> (random)] [--ene-min <min energy> (0.0)] [--ene-max <max energy> (1.0)] [--PE-steps <steps of PE evolution> (10)] [--thermalization <steps> (100)] [--record-reverse] [--walltime (0)]\n", argv[0]);
        exit(1);
    }

    parse_arguments(args, argc, argv);

    beta = args.beta;
    g_beta = args.g_beta; // defined as extern in system.cuh
    thermalization = args.thermalization;
    qms::metro_steps = (uint)args.metro_steps;
    qms::reset_each = (uint)args.reset_each;
    qms::state_qbits = (uint)args.state_qbits;
    qms::ene_qbits = (uint)args.ene_qbits;
    string outfilename(args.outfile);
    qms::max_reverse_attempts = (uint)args.max_reverse_attempts;
    qms::n_phase_estimation = args.pe_steps;
    qms::record_reverse= args.record_reverse;
    qms::iseed = args.seed;
    if(qms::iseed>0)
        qms::rangen.set_seed(qms::iseed);
    
    qms::iseed = qms::rangen.get_seed();

    qms::nqubits = qms::state_qbits + 2*qms::ene_qbits + 1;
    qms::Dim = (1U << qms::nqubits);
    qms::ene_levels = (1U << qms::ene_qbits);
    qms::state_levels = (1U << qms::state_qbits);

    qms::t_PE_shift = args.ene_min;
    qms::t_PE_factor = (qms::ene_levels-1)/(double)(qms::ene_levels*(args.ene_max-args.ene_min)); 
    qms::t_phase_estimation = qms::t_PE_factor*8.*atan(1.0); // 2*pi*t_PE_factor

    suqa::threads = NUM_THREADS;
    suqa::blocks = (qms::Dim+suqa::threads-1)/suqa::threads;
    if(suqa::blocks>MAXBLOCKS) suqa::blocks=MAXBLOCKS;

    
    // Banner
    suqa::print_banner();
    cout<<"arguments:\n"<<args<<endl;

    auto t_start = std::chrono::high_resolution_clock::now();

    // Initialization of utilities
    suqa::setup(qms::Dim);
    qms::setup(beta);

    // Initialization:
    // known eigenstate of the system (see src/system.cu)
    
    allocate_state(qms::gState, qms::Dim);
    init_state(qms::gState,qms::Dim);


    //TODO: make it an args option?
    uint perc_mstep = (qms::metro_steps+19)/20; // batched saves
    
    uint count_accepted = 0U;
//    if(!file_exists(outfilename.c_str())){
//        FILE * fil = fopen(outfilename.c_str(), "w");
//        fprintf(fil, "# E A\n");
//        fclose(fil);
//    }

    //XXX: systematic test
    int iiii=0;
    double rho_proj[8][8][2];
    
    string rhomat_fname="rho_mat_qms_b"+to_string(beta)+"_rt_"+to_string(qms::reset_each)+".txt";

    if( access( rhomat_fname.c_str(), F_OK ) != -1 ){
        printf("%s exists\n",rhomat_fname.c_str());
        FILE * fil = fopen(rhomat_fname.c_str(), "r");
        if(fscanf(fil,"%d\n",&iiii)!=1){ printf("Wrong reading!\n"); exit(1);}
        printf("Restarting from idx: %d\n",iiii);
        for(int i=0;i<8;++i){
            for(int j=0;j<8;++j) for(int k=0;k<2;++k){
                if(fscanf(fil,"%lg ",&rho_proj[i][j][k])!=1){ printf("Wrong reading!\n"); exit(1);}
            }
            if(fscanf(fil,"\n") !=0){ printf("Wrong reading!\n"); exit(1);}
        }
        fclose(fil);
        printf("Loading previous rho matrix");
    }else{
        for(uint i=0; i<8; ++i)  for(uint j=0; j<8; ++j) for(uint k=0; k<2; ++k) rho_proj[i][j][k]=0.0;
    }

    double TrDist_discrepancy, TrDist_discrepancy_prev=100000.;
    double Energy_discrepancy, Energy_discrepancy_prev=100000.;
    double Aoper_discrepancy, Aoper_discrepancy_prev=100000.;

    int TrDist_ctr=0;
    bool TrDist_first10done = false;
    std::vector<double> TrDist_recents(10);

    // partition function precomputation
    double Z=0.0;
    double E_sng_exact=0.0,E_sqr_exact;
    for(uint i=0; i<8; ++i){
        Z+=exp(-beta*eev[i]);
    }
    for(uint i=0;i<8;++i){
//        rho_diff_re[cci]=rho_proj[i][j][0]*rho_proj[i][j][0]/(sampling*sampling);
        E_sng_exact+=eev[i]*exp(-beta*eev[i])/Z;
        E_sqr_exact+=eev[i]*eev[i]*exp(-beta*eev[i])/Z;
    }

    auto t_prev = std::chrono::high_resolution_clock::now();

    bool take_measure;
    uint s0 = 0U;
    for(uint s = 0U; s < qms::metro_steps; ++s){
        DEBUG_CALL(cout<<"metro step: "<<s<<endl);
        take_measure = (s>s0+(uint)thermalization and (s-s0)%qms::reset_each ==0U);

        double tmp_rho[8][8][2];
        for(uint i=0; i<8; ++i)  for(uint j=0; j<8; ++j) for(uint k=0; k<2; ++k) tmp_rho[i][j][k]=0.0;

        int ret = qms::metro_step(take_measure, tmp_rho);

        // check conditions of measurement
        if(ret<0){ // failed rethermalization, reinitialize state
            init_state(qms::gState, qms::Dim);
            //ensure new rethermalization
            s0 = s+1; 
        }

        if(take_measure and (ret==2 or ret==4)){
            // measure rho as weighted average of eigenstates projectors
//            DUMP_STATE(qms::gState);

            iiii++;

//            printf("rho:\n");
            for(uint i=0;i<8;++i){
                for(uint j=0;j<8;++j){
                    rho_proj[i][j][0]+=tmp_rho[i][j][0];
                    rho_proj[i][j][1]+=tmp_rho[i][j][1];
//                    printf("(%.2lg %.2lg) ",rho_proj[i][j][0]/iiii,rho_proj[i][j][1]/iiii); 
                }
//                printf("\n");
            }

            auto t_tmp = std::chrono::high_resolution_clock::now();
            double secs_aft = (1./1000.)*std::chrono::duration<double, std::milli>(t_tmp-t_prev).count();
//            if(iiii%100==0){
            if(secs_aft>2.0){
                t_prev=t_tmp;
                double rho_diff_re[64]; //,rho_A_prod[8][8];
                double E_sng=0.0, E_sqr=0.0;
                double E_isolated=0.0;
                uint cci=0;
                for(uint i=0;i<8;++i) for(uint j=0;j<8;++j){
                //        rho_diff_re[cci]=rho_proj[i][j][0]*rho_proj[i][j][0]/(sampling*sampling);
                    rho_diff_re[cci]=rho_proj[i][j][0]/iiii;
                    //rho_A_prod[i][j]=0.0;

                    if(i==j){
                        E_isolated+=eev[i]*tmp_rho[i][j][0];
                        E_sng+=eev[i]*rho_diff_re[cci];
                        E_sqr+=eev[i]*eev[i]*rho_diff_re[cci];
                        
                        rho_diff_re[cci]-=exp(-beta*eev[i])/Z;
                    }
                //        printf("%.4lg %.4lg %.4lg\n",rho_proj[i][j][0]*rho_proj[i][j][0]/sampling,exp(-qsa::beta*eev[i])/Z, rho_diff_re[cci]);
                    cci++;
                }

//                printf("E_iso: %.6lg, E_meas: %.6lg\n",E_isolated,qms::E_measures.back());

                double A_sng=0.0, A_sqr=0.0;
                double A_sng_exact=0.0;
                for(int coeff=0;coeff<16;++coeff){
                    int i = Amat_coo[coeff][0];        
                    int k = Amat_coo[coeff][1];        
                    double a_ik= Amat_val[coeff][0]; // A_ik

                    for(int n=0; n<8; ++n){
                        for(int m=0; m<8; ++m){
                            A_sng+=ees[n][k]*rho_proj[n][m][0]/iiii*ees[m][i]*a_ik;
//                            A_sqr+=ees[n][k]*rho_proj[n][m][0]/iiii*ees[m][i]*a_ik*a_ik;
                        }
                        A_sng_exact+=ees[n][k]*exp(-beta*eev[n])/Z*ees[n][i]*a_ik;
                    }

                    for(int coeff2=0;coeff2<16;++coeff2){
                        int p = Amat_coo[coeff2][0];        
                        int l = Amat_coo[coeff2][1];        
                        if(p!=k) continue;
                        double a_kl= Amat_val[coeff][0]; // A_kl

                        for(int n=0; n<8; ++n) for(int m=0; m<8; ++m){
                                A_sqr+=ees[n][l]*rho_proj[n][m][0]/iiii*ees[m][i]*a_ik*a_kl;
                        }
                    }
                }

                printf("A_sng=%.8lg\tA_sqr=%.8lg; A_sng_exact=%.8lg\n",A_sng,A_sqr,A_sng_exact);
                double A_std = sqrt((A_sqr-A_sng*A_sng)/iiii);
                Aoper_discrepancy = (A_sng-A_sng_exact)/A_std;


                double E_std = sqrt((E_sqr-E_sng*E_sng)/iiii);
                double E_std_exact = sqrt((E_sqr_exact-E_sng_exact*E_sng_exact)/iiii);
                printf("E_exact: %.12lg, E: %.12lg, dE: %.12lg; rel_discr: %.12lg\n",E_sng_exact,E_sng,E_std_exact,abs(E_sng-E_sng_exact)/E_std_exact);

                Energy_discrepancy = (E_sng-E_sng_exact)/E_std;

                // eigensolver
                lapack_int lwork=3*8-1;
                double work[lwork];
                lapack_int n=8;
                lapack_int info;
                double dist_eigs[n];
                LAPACK_dsyev("N", "U", &n, rho_diff_re, &n, dist_eigs,work,&lwork,&info);

                if(info){
                    printf("Error in eigenvalue routine\n");
                }else{
                    double tr_dist=0.0;
//                    printf("dist eigs:\n");
                    for(uint i=0; i<8; ++i){
//                        printf("%.10lg\n",dist_eigs[i]);
                        tr_dist += abs(dist_eigs[i]);
                    }
                    tr_dist*=0.5;

//                    printf("Trace distance: %.12lg\n",tr_dist);
                    TrDist_discrepancy = tr_dist;
                }

                TrDist_recents[TrDist_ctr]=TrDist_discrepancy;
                TrDist_ctr=(TrDist_ctr+1)%10;
                if(!TrDist_first10done and TrDist_ctr==0) TrDist_first10done = true;

                double TrDist_ave=0.0, TrDist_fluct=0.0;
                if(TrDist_first10done){
                    for(int trdi=0;trdi<10;++trdi){
                        TrDist_ave+=TrDist_recents[trdi];
                        TrDist_fluct+=TrDist_recents[trdi]*TrDist_recents[trdi];
                    }
                    TrDist_ave/=10.0;
                    TrDist_fluct-=10.0*TrDist_ave*TrDist_ave;
                    TrDist_fluct*=10./9.;
                    TrDist_fluct=sqrt(TrDist_fluct);
                }


                double thr_discr=3.0;

                printf("%.8lg+-%.8lg (%.8lg);\t%.8lg+-%.8lg (%.8lg)|%.8lg\t%.8lg\t%.8lg\t%.8lg\n",E_sng,E_std,E_sng_exact,A_sng,A_std,A_sng_exact,TrDist_ave,TrDist_fluct,Energy_discrepancy,Aoper_discrepancy);

                if(    (abs(1.0-TrDist_discrepancy_prev/TrDist_discrepancy)<1e-3)
                   and (abs(Energy_discrepancy)>thr_discr or abs(1.0-Energy_discrepancy_prev/Energy_discrepancy)<1e-5)
                   and (abs(Aoper_discrepancy)>thr_discr or abs(1.0-Aoper_discrepancy_prev/Aoper_discrepancy)<1e-5)){
                    printf("Discrepancies converged\n");

                    if( access( outfilename.c_str(), F_OK ) == -1 ){
                        FILE * fil = fopen(outfilename.c_str(), "w");
                        //        fprintf(fil, "# E%s\n",(Xmatstem!="")?" A":"");
                        fprintf(fil, "#beta retherm E_mean E_err E_exact A_mean A_err A_exact TrDist TrDist_err\n");
                        fclose(fil);
                    }
                    FILE * fil = fopen(outfilename.c_str(), "a");
                    fprintf(fil, "%.8lg\t%2d\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\n", 
                            beta, qms::reset_each, E_sng, E_std, E_sng_exact, A_sng, A_std, A_sng_exact, TrDist_ave, TrDist_fluct);
                    fclose(fil);

                    break;
                }

                TrDist_discrepancy_prev = TrDist_discrepancy;
                Energy_discrepancy_prev = Energy_discrepancy;
                Aoper_discrepancy_prev = Aoper_discrepancy;

            }

        }

        if(ret==1 or ret==2){
            count_accepted++;
        }
        if(s%perc_mstep==0){
            cout<<"iteration: "<<s<<"/"<<qms::metro_steps<<endl;
//            save_measures(outfilename);
        }

        auto t_mid = std::chrono::high_resolution_clock::now();
        double secs_passed = (1./1000.)*std::chrono::duration<double, std::milli>(t_mid-t_start).count();
        if((args.walltime>0 and secs_passed>args.walltime) or access( "stop", F_OK ) != -1 ){
            remove("stop");
            // save rho and stop
            
            FILE * outrho = fopen(rhomat_fname.c_str(),"w");

            fprintf(outrho,"%d\n",iiii);
            for(int i=0;i<8;++i){
                for(int j=0;j<8;++j) for(int k=0;k<2;++k){
                    fprintf(outrho,"%.12lg ",rho_proj[i][j][k]);
                }
                fprintf(outrho,"\n");
            }
            fclose(outrho);

            printf("Closing due to walltime limit reached\n");
            break;
        }
    }

    cout<<endl;
    printf("\n\tacceptance: %3.2lg%%\n",(count_accepted/static_cast<double>(qms::metro_steps))*100.0);


    deallocate_state(qms::gState);
    qms::clear();
    suqa::clear();

    cout<<"\nall fine :)\n"<<endl;



    if(qms::record_reverse){
        FILE * fil_rev = fopen((outfilename+"_revcounts").c_str(), "w");

        for(uint i = 0; i < qms::reverse_counters.size(); ++i){
            fprintf(fil_rev, "%d %d\n", i, static_cast<int>(qms::reverse_counters[i]));
        }
        fclose(fil_rev);
    }

    cout<<"\n\tSuqa!\n"<<endl;

    auto t_end = std::chrono::high_resolution_clock::now();
    double secs_passed = (1./1000.)*std::chrono::duration<double, std::milli>(t_end-t_start).count();
	cout<<"All [DONE] in "<<secs_passed<<" seconds"<<endl;

    return 0;
}
