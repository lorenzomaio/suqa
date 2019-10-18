#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <cstring>
#include <stdio.h>
#include <bits/stdc++.h>
#include <cmath>
#include <cassert>
#include "include/Rand.hpp"
#include <chrono>

#ifndef NDEBUG
    #define DEBUG_CALL(x) x
#else
    #define DEBUG_CALL(x)
#endif

using namespace std;

void print_banner(){
    printf("\n"
"                                          \n" 
"    ███████╗██╗   ██╗ ██████╗  █████╗     \n" 
"    ██╔════╝██║   ██║██╔═══██╗██╔══██╗    \n" 
"    ███████╗██║   ██║██║   ██║███████║    \n" 
"    ╚════██║██║   ██║██║▄▄ ██║██╔══██║    \n" 
"    ███████║╚██████╔╝╚██████╔╝██║  ██║    \n" 
"    ╚══════╝ ╚═════╝  ╚══▀▀═╝ ╚═╝  ╚═╝    \n" 
"                                          \n" 
"\nSimulator for Universal Quantum Algorithms\n"); 
}


typedef complex<double> Complex;
const Complex iu(0, 1);

/* Hamiltonian
 *
 * H = {{1, 0, 0}, {0, 1, 1}, {0, 1, 1}}
 *
 */


// simulation parameters
double beta;
double f1;
double f2;
const uint nqubits = 7;
const uint Dim = (uint)pow(2.0, nqubits); // simulation hyperparameters
uint max_reverse_attempts;
uint metro_steps;
uint reset_each;
unsigned long long iseed = 0ULL;
double t_phase_estimation;
int n_phase_estimation;

uint gCi;
uint c_acc = 0;


// Global state of the system.
// Ordering (less to most significant)
// psi[0], psi[1], E_old[0], E_old[1], E_new[0], E_new[1], acc //, qaux[0]
vector<Complex> gState(Dim,0.0);

//vector<double> energy_measures;
vector<double> X_measures;
vector<double> E_measures;

// Operator X parameter
const double phi = (1.+sqrt(5.))/2.;
const double mphi_inv = -1./phi;
const double Sa = phi/sqrt(2.+phi);
const double Sb = 1/sqrt(2.+phi);
const double S_10=Sa, S_12=Sb, S_20=-Sb, S_22=Sa;
const double twosqinv = 1./sqrt(2.);


// constants
const Complex rphase_m[3] = {exp((2*M_PI)*iu), exp((2*M_PI/2.)*iu), exp((2*M_PI/4.)*iu)};

// Utilities

pcg rangen;

// bit masks
enum bm_idxs {  bm_psi0, 
                bm_psi1,
                bm_E_old0,
                bm_E_old1,
                bm_E_new0,
                bm_E_new1,
                bm_acc};


std::ostream& operator<<(std::ostream& s, const Complex& c){
    s<<"("<<real(c)<<", "<<imag(c)<<")";
    return s;
}

template<class T>
void print(vector<T> v){
    for(const auto& el : v)
        cout<<el<<" ";
    cout<<endl;
}

void sparse_print(vector<Complex> v){
    for(uint i=0; i<v.size(); ++i){
        if(norm(v[i])>1e-8)
            cout<<"i="<<i<<" -> "<<v[i]<<"; ";
    }
    cout<<endl;
}

double vnorm(const vector<Complex>& v){
    double ret = 0.0;
    for(const auto& el : v)
        ret += norm(el);
    
    return sqrt(ret);
}

void vnormalize(vector<Complex>& v){
    double vec_norm = vnorm(v);
    for(auto& el : v)
        el/=vec_norm;
}

template<class T>
void apply_2x2mat(T& x1, T& x2, const Complex& m11, const Complex& m12, const Complex& m21, const Complex& m22){

            T x1_next = m11 * x1 + m12 * x2;
            T x2_next = m21 * x1 + m22 * x2;
            x1 = x1_next;
            x2 = x2_next;
}

void qi_reset(vector<Complex>& state, const uint& q){
    for(uint i = 0U; i < state.size(); ++i){
        if((i >> q) & 1U){ // checks q-th digit in i
            uint j = i & ~(1U << q); // j has 0 on q-th digit
            state[j]+=state[i];
            state[i]= {0.0, 0.0};
        }
    }
    vnormalize(state);
}  

void qi_reset(vector<Complex>& state, const vector<uint>& qs){
    uint mask=0U;
    for(const auto& q : qs)
        mask |= 1U << q;

    for(uint i = 0U; i < state.size(); ++i){
        if((i & mask) != 0U){ // checks q-th digit in i
            state[i & ~mask]+=state[i];
            state[i]= 0.0;
        }
    }
    vnormalize(state);
}  

void qi_x(vector<Complex>& state, const uint& q){
    for(uint i = 0U; i < state.size(); ++i){
        if((i >> q) & 1U){ // checks q-th digit in i
            uint j = i & ~(1U << q); // j has 0 on q-th digit
            std::swap(state[i],state[j]);
        }
    }
}  

void qi_x(vector<Complex>& state, const vector<uint>& qs){
    for(const auto& q : qs)
        qi_x(state, q);
}  

void qi_h(vector<Complex>& state, const uint& q){
	for(uint i_0 = 0U; i_0 < state.size(); ++i_0){
        if((i_0 & (1U << q)) == 0U){
            uint i_1 = i_0 | (1U << q);
            Complex a_0 = state[i_0];
            Complex a_1 = state[i_1];
            state[i_0] = twosqinv*(a_0+a_1);
            state[i_1] = twosqinv*(a_0-a_1);
        }
    }
}  

void qi_h(vector<Complex>& state, const vector<uint>& qs){
    for(const auto& q : qs)
        qi_h(state, q);
}  


void qi_cx(vector<Complex>& state, const uint& q_control, const uint& q_target){
    for(uint i = 0U; i < state.size(); ++i){
        // for the swap, not only q_target:1 but also q_control:1
        if(((i >> q_control) & 1U) && ((i >> q_target) & 1U)){
            uint j = i & ~(1U << q_target);
            std::swap(state[i],state[j]);
        }
    }
}  
  

void qi_cx(vector<Complex>& state, const uint& q_control, const uint& q_mask, const uint& q_target){
    uint mask_qs = 1U << q_target;
    uint mask = mask_qs | (1U << q_control);
    if(q_mask) mask_qs |= (1U << q_control);
    for(uint i = 0U; i < state.size(); ++i){
        if((i & mask) == mask_qs){
            uint j = i & ~(1U << q_target);
            std::swap(state[i],state[j]);
        }
    }
}  

void qi_mcx(vector<Complex>& state, const vector<uint>& q_controls, const uint& q_target){
    uint mask = 1U << q_target;
    for(const auto& q : q_controls)
        mask |= 1U << q;

    for(uint i = 0U; i < state.size(); ++i){
        if((i & mask) == mask){
            uint j = i & ~(1U << q_target);
            std::swap(state[i],state[j]);
        }
    }
}  


void qi_mcx(vector<Complex>& state, const vector<uint>& q_controls, const vector<uint>& q_mask, const uint& q_target){
    uint mask = 1U << q_target;
    for(const auto& q : q_controls)
        mask |= 1U << q;
    uint mask_qs = 1U << q_target;
    for(uint k = 0U; k < q_controls.size(); ++k){
        if(q_mask[k]) mask_qs |= 1U << q_controls[k];
    }
    for(uint i = 0U; i < state.size(); ++i){
        if((i & mask) == mask_qs){
            uint j = i & ~(1U << q_target);
            std::swap(state[i],state[j]);
        }
    }
}  
void qi_swap(vector<Complex>& state, const uint& q1, const uint& q2){
        // swap gate: 00->00, 01->10, 10->01, 11->11
        // equivalent to cx(q1,q2)->cx(q2,q1)->cx(q1,q2)
        uint mask_q = (1U << q1);
        uint mask = mask_q | (1U << q2);
        for(uint i = 0U; i < state.size(); ++i){
            if((i & mask) == mask_q){
                uint j = (i & ~(1U << q1)) | (1U << q2);
                std::swap(state[i],state[j]);
            }
        }
}

// Simulation procedures

void check_unused(){
    uint mask1 = 3U;
    uint mask2 = 12U;
    for(uint i = 0U; i < gState.size(); ++i){
        if( (i & mask1) == mask1){
            assert(norm(gState[i])<1e-8);
        }
        if( (i & mask2) == mask2){
            assert(norm(gState[i])<1e-8);
        }
    } 
}

void reset_non_state_qbits(){
    DEBUG_CALL(cout<<"\n\nBefore reset"<<endl);
    DEBUG_CALL(sparse_print(gState));
    qi_reset(gState, {bm_E_old0, bm_E_old1, bm_E_new0, bm_E_new1, bm_acc});
    DEBUG_CALL(cout<<"\n\nAfter reset"<<endl);
    DEBUG_CALL(sparse_print(gState));
}


void measure_qbit(vector<Complex>& state, const uint& q, uint& c){
    double prob1 = 0.0;

    for(uint i = 0U; i < state.size(); ++i){
        if((i >> q) & 1U){
            prob1+=norm(state[i]); 
        }
    }
    c = (uint)(rangen.doub() < prob1); // prob1=1 -> c = 1 surely
    
    if(c){ // set to 0 coeffs with bm_acc 0
        for(uint i = 0U; i < state.size(); ++i){
            if(((i >> q) & 1U) == 0U)
                state[i] = {0.0, 0.0};        
        }
    }else{ // set to 0 coeffs with bm_acc 1
        for(uint i = 0U; i < state.size(); ++i){
            if(((i >> q) & 1U) == 1U)
                state[i] = {0.0, 0.0};        
        }
    }
    vnormalize(state);
}

//TODO: can be optimized for multiple qbits measures?
void measure_qbits(vector<Complex>& state, const vector<uint>& qs, vector<uint>& cs){
    for(uint k = 0U; k < qs.size(); ++k)
        measure_qbit(state, qs[k], cs[k]);
}

void qi_crm(vector<Complex>& state, const uint& q_control, const uint& q_target, const int& m){
    for(uint i = 0U; i < state.size(); ++i){
        // for the swap, not only q_target:1 but also q_control:1
        if(((i >> q_control) & 1U) && ((i >> q_target) & 1U)){
            state[i] *= (m>0) ? rphase_m[m] : conj(rphase_m[-m]);
        }
    }
}

void qi_cu_on2(vector<Complex>& state, const double& dt, const uint& q_control, const vector<uint>& qstate){
    uint cmask = (1U << q_control);
	uint mask = cmask; // (1U << qstate[0]) | (1U << qstate[0])
    for(const auto& qs : qstate){
        mask |= (1U << qs);
    }

	for(uint i_0 = 0U; i_0 < state.size(); ++i_0){
        if((i_0 & mask) == cmask){
      
            uint i_1 = i_0 | (1U << qstate[0]);
            uint i_2 = i_0 | (1U << qstate[1]);
            uint i_3 = i_1 | i_2;

            Complex a_0 = state[i_0];
            Complex a_1 = state[i_1];
            Complex a_2 = state[i_2];
            Complex a_3 = state[i_3];
            

            state[i_0] = exp(-dt*iu)*a_0;
            state[i_1] = exp(-dt*iu)*(cos(dt)*a_1 -sin(dt)*iu*a_2);
            state[i_2] = exp(-dt*iu)*(-sin(dt)*iu*a_1 + cos(dt)*a_2);
            state[i_3] = exp(-dt*iu)*a_3;
        }
    }

}

void qi_qft(vector<Complex>& state, const vector<uint>& qact){
    if(qact.size()!=2)
        throw "ERROR: qft(inverse) not implemented for nqubits != 2";

    qi_h(state, qact[1]);
    qi_crm(state, qact[0], qact[1], -2);
    qi_h(state, qact[0]);

}

void qi_qft_inverse(vector<Complex>& state, const vector<uint>& qact){
    if(qact.size()!=2)
        throw "ERROR: qft(inverse) not implemented for nqubits != 2";

    qi_h(state, qact[0]);
    qi_crm(state, qact[0], qact[1], 2);
    qi_h(state, qact[1]);

}

void apply_phase_estimation(vector<Complex>& state, const vector<uint>& q_state, const vector<uint>& q_target, const double& t, const uint& n){
    DEBUG_CALL(cout<<"apply_phase_estimation()"<<endl);
    qi_h(state,q_target);
    DEBUG_CALL(cout<<"after qi_h(state,q_target)"<<endl);
    DEBUG_CALL(sparse_print(state));

    // apply CUs
    double dt = t/(double)n;

    for(int trg = q_target.size() - 1; trg > -1; --trg){
        for(uint ti = 0; ti < n; ++ti){
            for(uint itrs = 0; itrs < q_target.size()-trg; ++itrs){
                qi_cu_on2(state, dt, q_target[trg], q_state);
            }
        }
    }
    DEBUG_CALL(cout<<"\nafter evolutions"<<endl);
    DEBUG_CALL(sparse_print(state));
    
    // apply QFT^{-1}
    qi_qft_inverse(state, q_target); 

}

void apply_phase_estimation_inverse(vector<Complex>& state, const vector<uint>& q_state, const vector<uint>& q_target, const double& t, const uint& n){
    DEBUG_CALL(cout<<"apply_phase_estimation_inverse()"<<endl);

    // apply QFT
    qi_qft(state, q_target); 


    // apply CUs
    double dt = t/(double)n;

    for(uint trg = 0; trg < q_target.size(); ++trg){
        for(uint ti = 0; ti < n; ++ti){
            for(uint itrs = 0; itrs < q_target.size()-trg; ++itrs){
                qi_cu_on2(state, -dt, q_target[trg], q_state);
            }
        }
    }
    
    qi_h(state,q_target);

}

// void apply_Phi_old(){
//    // quantum phase estimation (here trivial)
//     DEBUG_CALL(cout<<"\nApply Phi old\n"<<endl);
// 	uint mask = 3U;
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 - twosqinv*a_2;
//             gState[i_2] = twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
// 
//     qi_mcx(gState, {bm_psi0,bm_psi1}, {0,0},bm_E_old0);
//     qi_cx(gState, bm_psi1, bm_E_old1);
// 
//     DEBUG_CALL(cout<<"\napply Phi_diag old:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 + twosqinv*a_2;
//             gState[i_2] = -twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S_dag old:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// }
// 
// void apply_Phi_old_inverse(){
//    // quantum phase estimation (here trivial)
//     DEBUG_CALL(cout<<"\nApply Phi old inverse\n"<<endl);
// 	uint mask = 3U;
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 - twosqinv*a_2;
//             gState[i_2] = twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
//     qi_mcx(gState, {bm_psi0,bm_psi1}, {0,0},bm_E_old0);
//     qi_cx(gState, bm_psi1, bm_E_old1);
// 
//     DEBUG_CALL(cout<<"\napply Phi_diag old inverse:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 + twosqinv*a_2;
//             gState[i_2] = -twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S_dag old inverse:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// }
// 
// void apply_Phi(){
//    // quantum phase estimation (here trivial)
//     DEBUG_CALL(cout<<"\nApply Phi\n"<<endl);
// 	uint mask = 3U;
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 - twosqinv*a_2;
//             gState[i_2] = twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
//     qi_mcx(gState, {bm_psi0,bm_psi1}, {0,0},bm_E_new0);
//     qi_cx(gState, bm_psi1, bm_E_new1);
//     
//     DEBUG_CALL(cout<<"\napply Phi_diag:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 + twosqinv*a_2;
//             gState[i_2] = -twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S_dag:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// }
// 
// void apply_Phi_inverse(){
//    // quantum phase estimation (here trivial)
//     DEBUG_CALL(cout<<"\nApply Phi inverse\n"<<endl);
// 	uint mask = 3U;
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 - twosqinv*a_2;
//             gState[i_2] = twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
//     qi_mcx(gState, {bm_psi0,bm_psi1}, {0,0},bm_E_new0);
//     qi_cx(gState, bm_psi1, bm_E_new1);
// 
//     DEBUG_CALL(cout<<"\napply Phi_diag inverse:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// 
// 	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
//         if((i_0 & mask) == 0U){
//             uint i_1 = i_0 | 1U;
//             uint i_2 = i_0 | 2U;
// 
//             Complex a_0 = gState[i_0];
//             Complex a_1 = gState[i_1];
//             Complex a_2 = gState[i_2];
//             
//             gState[i_0] = a_0;
//             gState[i_1] = twosqinv*a_1 + twosqinv*a_2;
//             gState[i_2] = -twosqinv*a_1 + twosqinv*a_2;
//         }
//     }
//     DEBUG_CALL(cout<<"\napply S_dag:\n"<<endl);
//     DEBUG_CALL(sparse_print(gState));
// }



void apply_Phi_old(){

    apply_phase_estimation(gState, {bm_psi0, bm_psi1}, {bm_E_old0, bm_E_old1}, t_phase_estimation, n_phase_estimation);

}

void apply_Phi_old_inverse(){

    apply_phase_estimation_inverse(gState, {bm_psi0, bm_psi1}, {bm_E_old0, bm_E_old1}, t_phase_estimation, n_phase_estimation);

}

void apply_Phi(){

    apply_phase_estimation(gState, {bm_psi0, bm_psi1}, {bm_E_new0, bm_E_new1}, t_phase_estimation, n_phase_estimation);

}

void apply_Phi_inverse(){

    apply_phase_estimation_inverse(gState, {bm_psi0, bm_psi1}, {bm_E_new0, bm_E_new1}, t_phase_estimation, n_phase_estimation);

}


uint draw_C(){
    if (rangen.doub()<0.5)
        return 0U;
    return 1U;
}

void apply_C(const uint &Ci){
    if(Ci==0U){
        qi_cx(gState,bm_psi1, 0, bm_psi0);
    }else if(Ci==1U){
        qi_swap(gState,bm_psi1,bm_psi0);
    }else{
        throw "Error!";
    }
}

void apply_C_inverse(const uint &Ci){
    apply_C(Ci);
}

void apply_W(){
    DEBUG_CALL(cout<<"\n\nApply W"<<endl);
    //(1U <<bm_E_old0) | (1U <<bm_E_old1) |(1U <<bm_E_new0) |(1U <<bm_E_new1) |(1U <<bm_acc); 
    uint mask = 124U;
    // Ei = 0, Ek = 1
    //(1U <<bm_E_new0) |(1U <<bm_acc);
    uint case1a = 80U;
    // Ei = 1, Ek = 2
    //(1U <<bm_E_old0) |(1U <<bm_E_new1) |(1U <<bm_acc);
    uint case1b = 100U;
    // Ei = 0, Ek = 2
    //(1U <<bm_E_new1) |(1U <<bm_acc);
    uint case2 = 96U;
    for(uint i = 0U; i < gState.size(); ++i){
        if(((i & mask) == case1a) || ((i & mask) == case1b)){
            uint j = i & ~(1U << bm_acc);
            
            DEBUG_CALL(if(norm(gState[i])+norm(gState[j])>1e-8) cout<<"case1: gState["<<i<<"] = "<<gState[i]<<", gState["<<j<<"] = "<<gState[j]<<endl);
            apply_2x2mat(gState[j], gState[i], sqrt(1.-f1), sqrt(f1), sqrt(f1), -sqrt(1.-f1));
            DEBUG_CALL(if(norm(gState[i])+norm(gState[j])>1e-8) cout<<"after: gState["<<i<<"] = "<<gState[i]<<", gState["<<j<<"] = "<<gState[j]<<endl);
        }else if((i & mask) == case2){
            uint j = i & ~(1U << bm_acc);

            DEBUG_CALL(if(norm(gState[i])+norm(gState[j])>1e-8) cout<<"case2: gState["<<i<<"] = "<<gState[i]<<", gState["<<j<<"] = "<<gState[j]<<endl);
            apply_2x2mat(gState[j], gState[i], sqrt(1.-f2), sqrt(f2), sqrt(f2), -sqrt(1.-f2));
            DEBUG_CALL(if(norm(gState[i])+norm(gState[j])>1e-8) cout<<"after: gState["<<i<<"] = "<<gState[i]<<", gState["<<j<<"] = "<<gState[j]<<endl);
        }else if((i >> bm_acc) & 1U){
            uint j = i & ~(1U << bm_acc);

            DEBUG_CALL(if(norm(gState[i])+norm(gState[j])>1e-8) cout<<"case3: gState["<<i<<"] = "<<gState[i]<<", gState["<<j<<"] = "<<gState[j]<<endl);
            std::swap(gState[i],gState[j]);
            DEBUG_CALL(if(norm(gState[i])+norm(gState[j])>1e-8) cout<<"after: gState["<<i<<"] = "<<gState[i]<<", gState["<<j<<"] = "<<gState[j]<<endl);
        }
    }
}

void apply_W_inverse(){
    apply_W();
}

void apply_U(){
    DEBUG_CALL(cout<<"\n\nApply U"<<endl);
    apply_C(gCi);
    DEBUG_CALL(cout<<"\n\nAfter apply C = "<<gCi<<endl);
    DEBUG_CALL(sparse_print(gState));
    apply_Phi();
    DEBUG_CALL(cout<<"\n\nAfter second phase estimation"<<endl);
    DEBUG_CALL(sparse_print(gState));
    apply_W();
    DEBUG_CALL(cout<<"\n\nAfter apply W"<<endl);
    DEBUG_CALL(sparse_print(gState));
}

void apply_U_inverse(){
    apply_W_inverse();
    DEBUG_CALL(cout<<"\n\nAfter apply W inverse"<<endl);
    DEBUG_CALL(sparse_print(gState));
    apply_Phi_inverse();
    DEBUG_CALL(cout<<"\n\nAfter inverse second phase estimation"<<endl);
    DEBUG_CALL(sparse_print(gState));
    apply_C_inverse(gCi);
    DEBUG_CALL(cout<<"\n\nAfter apply C inverse = "<<gCi<<endl);
    DEBUG_CALL(sparse_print(gState));
}

double measure_X(){
	uint mask = 3U;
	vector<uint> classics(2);
	for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
        if((i_0 & mask) == 0U){
            uint i_1 = i_0 | 1U;
            uint i_2 = i_0 | 2U;

            Complex a_0 = gState[i_0];
            Complex a_1 = gState[i_1];
            Complex a_2 = gState[i_2];
            
            gState[i_0] = a_1;
            gState[i_1] = Sa*a_0 + Sb*a_2;
            gState[i_2] = -Sb*a_0 + Sa*a_2;
        }
    }
    measure_qbits(gState, {bm_psi0,bm_psi1}, classics);
    for(uint i_0 = 0U; i_0 < gState.size(); ++i_0){
        if((i_0 & mask) == 0U){
            uint i_1 = i_0 | 1U;
            uint i_2 = i_0 | 2U;

            Complex a_0 = gState[i_0];
            Complex a_1 = gState[i_1];
            Complex a_2 = gState[i_2];

            gState[i_0] = Sa*a_1 - Sb*a_2;
            gState[i_1] = a_0;
            gState[i_2] = Sb*a_1 + Sa*a_2;
        }
    }
    uint meas = classics[0] + 2*classics[1];
    switch(meas){
        case 0:
            return 0;
            break;
        case 1:
            return phi;
            break;
        case 2:
            return mphi_inv;
            break;
        default:
            throw "Error!";
    }
    return 0.0;
}

// double measure_X(){
// 	vector<uint> classics(2);
//     measure_qbits(gState, {bm_psi0,bm_psi1}, classics);
//     uint meas = classics[0] + 2*classics[1];
//     switch(meas){
//         case 0:
//             return 1.0;
//             break;
//         case 1:
//             return 2.0;
//             break;
//         case 2:
//             return 3.0;
//             break;
//         default:
//             throw "Error!";
//     }
//     return 0.0;
// }


void metro_step(uint s){
    DEBUG_CALL(cout<<"initial state"<<endl);
    DEBUG_CALL(sparse_print(gState));
    reset_non_state_qbits();
    DEBUG_CALL(check_unused());
    DEBUG_CALL(cout<<"state after reset"<<endl);
    DEBUG_CALL(sparse_print(gState));
    apply_Phi_old();
    DEBUG_CALL(check_unused());
    DEBUG_CALL(cout<<"\n\nAfter first phase estimation"<<endl);
    DEBUG_CALL(sparse_print(gState));

    gCi = draw_C();
    DEBUG_CALL(cout<<"\n\ndrawn C = "<<gCi<<endl);
    apply_U();
    DEBUG_CALL(check_unused());

    measure_qbit(gState, bm_acc, c_acc);

    if (c_acc == 1U){
        DEBUG_CALL(cout<<"accepted"<<endl);
        vector<uint> c_E_news(2,0), c_E_olds(2,0);
        measure_qbits(gState, {bm_E_new0, bm_E_new1}, c_E_news);
        DEBUG_CALL(double tmp_E=c_E_news[0]+2*c_E_news[1]);
        DEBUG_CALL(cout<<"  energy measure : "<<tmp_E<<endl); 
        apply_Phi_inverse();
//        E_measures.push_back(c_E_news[0]+2*c_E_news[1]);
        if(s>0U and s%reset_each ==0U){
            E_measures.push_back(c_E_news[0]+2*c_E_news[1]);
            qi_reset(gState, {bm_E_new0, bm_E_new1});
            X_measures.push_back(measure_X());
////            X_measures.push_back(0.0);
            DEBUG_CALL(cout<<"  X measure : "<<X_measures.back()<<endl); 
            DEBUG_CALL(cout<<"\n\nAfter X measure"<<endl);
            DEBUG_CALL(sparse_print(gState));
            DEBUG_CALL(cout<<"  X measure : "<<X_measures.back()<<endl); 
//            reset_non_state_qbits();
            qi_reset(gState, {bm_E_new0, bm_E_new1});
            apply_Phi();
            measure_qbits(gState, {bm_E_new0, bm_E_new1}, c_E_news);
            DEBUG_CALL(cout<<"\n\nAfter E recollapse"<<endl);
            DEBUG_CALL(sparse_print(gState));
            apply_Phi_inverse();
      }

        return;
    }
    //else

    DEBUG_CALL(cout<<"rejected; restoration cycle:"<<endl);
    apply_U_inverse();

    DEBUG_CALL(cout<<"\n\nBefore reverse attempts"<<endl);
    DEBUG_CALL(sparse_print(gState));
    uint iters = max_reverse_attempts;
    while(iters > 0){
        apply_Phi();
        double Eold_meas, Enew_meas;
        vector<uint> c_E_olds(2,0);
        vector<uint> c_E_news(2,0);
        measure_qbits(gState, {bm_E_old0, bm_E_old1}, c_E_olds);
        Eold_meas = c_E_olds[0]+2*c_E_olds[1];
        measure_qbits(gState, {bm_E_new0, bm_E_new1}, c_E_news);
        Enew_meas = c_E_news[0]+2*c_E_news[1];
        apply_Phi_inverse();
        
        if(Eold_meas == Enew_meas){
//            E_measures.push_back(Eold_meas);
            DEBUG_CALL(cout<<"  accepted restoration ("<<max_reverse_attempts-iters<<"/"<<max_reverse_attempts<<")"<<endl); 
            if(s>0U and s%reset_each == 0U){
                E_measures.push_back(Eold_meas);
                DEBUG_CALL(cout<<"  energy measure : "<<Eold_meas<<endl); 
                DEBUG_CALL(cout<<"\n\nBefore X measure"<<endl);
                DEBUG_CALL(sparse_print(gState));
                qi_reset(gState, {bm_E_new0, bm_E_new1});
                X_measures.push_back(measure_X());
////                X_measures.push_back(0.);
                DEBUG_CALL(cout<<"\n\nAfter X measure"<<endl);
                DEBUG_CALL(sparse_print(gState));
                DEBUG_CALL(cout<<"  X measure : "<<X_measures.back()<<endl); 
 ////               reset_non_state_qbits();
                qi_reset(gState, {bm_E_new0, bm_E_new1});
                apply_Phi();
                measure_qbits(gState, {bm_E_new0, bm_E_new1}, c_E_news);
                DEBUG_CALL(cout<<"\n\nAfter E recollapse"<<endl);
                DEBUG_CALL(sparse_print(gState));
                apply_Phi_inverse();
            }
            break;
        }
        //else
        DEBUG_CALL(cout<<"  rejected ("<<max_reverse_attempts-iters<<"/"<<max_reverse_attempts<<")"<<endl); 
        uint c_acc_trash;
        apply_U(); 
        measure_qbit(gState, bm_acc, c_acc_trash); 
        apply_U_inverse(); 

        iters--;
    }

    if (iters == 0){
        cout<<"not converged :("<<endl;
        exit(1);
    }
}


struct arg_list{
    double beta = 0.0;
    int metro_steps = 0;
    int reset_each = 0;
    string outfile = "";
    int max_reverse_attempts = 100;
    unsigned long long int seed = 0;
    double pe_time = 2.*atan(1.0);
    int pe_steps = 10; 
    

    friend ostream& operator<<(ostream& o, const arg_list& al);
};

ostream& operator<<(ostream& o, const arg_list& al){
    o<<"beta: "<<al.beta<<endl;
    o<<"metro steps: "<<al.metro_steps<<endl;
    o<<"reset each: "<<al.reset_each<<endl;
    o<<"max reverse attempts: "<<al.max_reverse_attempts<<endl;
    o<<"seed: "<<al.seed<<endl;
    o<<"out datafile: "<<al.outfile<<endl;
    o<<"time of PE evolution: "<<al.pe_time<<endl;
    o<<"steps of PE evolution: "<<al.pe_steps<<endl;
    return o;
}

arg_list args;

void parse_arguments(arg_list& args, int argc, char** argv){
    int fixed_args = 4;
    map<string,int> argmap;
    map<int,string> argmap_inv;
    char *end;
    int base_strtoull = 10;

    // fixed arguments
    args.beta = stod(argv[1],NULL);
    args.metro_steps = atoi(argv[2]);
    args.reset_each = atoi(argv[3]);
    args.outfile = argv[4];

    // floating arguments
    for(int i = fixed_args+1; i < argc; ++i){
        argmap[argv[i]]=i;
        argmap_inv[i]=argv[i];
    }
    int tmp_idx;

    // (int) max_reverse_attempts
    tmp_idx = argmap["--max-reverse"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--max-reverse' flag"; 
       
       args.max_reverse_attempts = atoi(argmap_inv[tmp_idx+1].c_str()); 
    }

    // (unsigned long long) seed
    tmp_idx = argmap["--seed"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--seed' flag"; 
       
       args.seed = strtoull(argmap_inv[tmp_idx+1].c_str(), &end, base_strtoull); 
    }

    // (double) pe_time
    tmp_idx = argmap["--PE-time"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--PE-time' flag"; 
       
       args.pe_time *= stod(argmap_inv[tmp_idx+1].c_str(), NULL); 
    }

    // (int) pe_steps
    tmp_idx = argmap["--PE-steps"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--PE-steps' flag"; 
       
       args.pe_steps = stod(argmap_inv[tmp_idx+1].c_str(), NULL); 
    }


    // argument checking
    if(args.beta <= 0.0){
        throw "ERROR: argument <beta> invalid";
    }

    if(args.metro_steps <= 0){
        throw "ERROR: argument <metro steps> invalid";
    }

    if(args.reset_each <=0){
        throw "ERROR: argument <reset each> non positive";
    }
    
    if(args.outfile == ""){
        throw "ERROR: argument <output file path> empty";
    }

    if(args.max_reverse_attempts <=0){
        throw "ERROR: argument <max reverse attempts> non positive";
    }

    if(args.pe_time <=0){
        throw "ERROR: argument <time of PE evolution> non positive";
    }

    if(args.pe_steps <=0){
        throw "ERROR: argument <steps of PE evolution> non positive";
    }
}


int main(int argc, char** argv){
    if(argc < 5){
        cout<<"arguments: <beta> <metro steps> <reset each> <output file path> [--max-reverse <max reverse attempts>=20] [--seed <seed>=random] [--PE-time <time of PE evolution>=pi/2] [--PE-steps <steps of PE evolution>=10]"<<endl;
        exit(1);
    }

    parse_arguments(args, argc, argv);

    beta = args.beta;
    metro_steps = (uint)args.metro_steps;
    reset_each = (uint)args.reset_each;
    string outfilename(args.outfile);
    max_reverse_attempts = (uint)args.max_reverse_attempts;
    t_phase_estimation = args.pe_time;
    n_phase_estimation = args.pe_steps;
    iseed = args.seed;
    if(iseed>0)
        rangen.set_seed(iseed);
    
    iseed = rangen.get_seed();

    f1 = exp(-beta);
    f2 = exp(-2.*beta);

    
    // Banner
    print_banner();
    cout<<"arguments:\n"<<args<<endl;

    auto t_start = std::chrono::high_resolution_clock::now();

    // Initialization:
    // known eigenstate of the system: psi=0, E_old = 0
    
    std::fill_n(gState.begin(), gState.size(), 0.0);
    gState[0] = 1.0; 
    for(uint s = 0U; s < metro_steps; ++s){
        metro_step(s);
    }

    cout<<"all fine :)\n"<<endl;

    FILE * fil = fopen(outfilename.c_str(), "w");

    fprintf(fil, "# it E X\n");

    for(uint ei = 0; ei < E_measures.size(); ++ei){
        fprintf(fil, "%d %.16lg %.16lg\n", ei, E_measures[ei], X_measures[ei]);
    }
//    for(uint ei = 0; ei < E_measures.size(); ++ei){
//        fprintf(fil, "%d %.16lg\n", ei, E_measures[ei]);
//    }
    fclose(fil);

    cout<<"\n\tSuqa!\n"<<endl;

    // test gates:
//    {
//        cout<<"TEST GATES"<<endl;
//        vector<Complex> test_state = {{0.4,-1.6},{1.2,0.7},{-0.1,0.6},{-1.3,0.4},{1.2,-1.3},{-1.2,1.7},{-3.1,0.6},{-0.3,0.2}};
//        vnormalize(test_state);
//        cout<<"initial state:"<<endl;
//        print(test_state);
//        cout<<"apply X to qbit 1 (most significant one)"<<endl;
//        qi_x(test_state, 1);
//        print(test_state);
//        cout<<"reapply X to qbit 1 (most significant one)"<<endl;
//        qi_x(test_state, 1);
//        print(test_state);
//        cout<<"apply CX controlled by qbit 0 to qbit 1"<<endl;
//        qi_cx(test_state, 0, 1);
//        print(test_state);
//        cout<<"apply CCX controlled by qbit 1 and 2 to qbit 0"<<endl;
//        qi_mcx(test_state, {1,2}, 0);
//        print(test_state);
//    }
//    { 
//        cout<<"\nTEST SIMULATION"<<endl;
//        vector<Complex> test_state = {{0.4,-1.6},{1.2,0.7},{-0.1,0.6},{-1.3,0.4},{1.2,-1.3},{-1.2,1.7},{-3.1,0.6},{-0.3,0.2}};
//        vnormalize(test_state);
//        gState = test_state;
//        cout<<"initial state:"<<endl;
//        sparse_print(gState);
//        for(uint jj=0; jj<3; ++jj){
//            cout<<"draw C (qubits 0 and 1 involved)"<<endl;
//            gCi = draw_C();
//            cout<<"drawn "<<gCi<<", apply it"<<endl;
//            apply_C(gCi);
//            sparse_print(gState);
//        }
//        cout<<"measure qubit 1"<<endl;
//        uint ctest;
//        measure_qbit(gState, 1U, ctest);
//        sparse_print(gState);
//    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double secs_passed = (1./1000.)*std::chrono::duration<double, std::milli>(t_end-t_start).count();
	cout<<"All [DONE] in "<<secs_passed<<" seconds"<<endl;

    return 0;
}
