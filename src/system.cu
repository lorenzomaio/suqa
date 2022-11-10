#include "system.cuh"

double g_beta;

int num_phys_levels;
std::vector<double> phys_spectrum;

void init_state(){
    suqa::init_state();
}


/* Quantum evolutor of the state */
void evolution(const double& t, const int& n){
    (void)n;

    std::vector<double> evophases(1U<<syst_qbits,0.);
    for(int ii=0; ii<num_phys_levels; ++ii){
        evophases[ii] = -phys_spectrum[ii]*t;
    }

    suqa::apply_phase_list(0,(uint)syst_qbits,evophases);
}



// qsa specifics
void qsa_init_state(){
    //TODO: implement
    throw std::runtime_error("ERROR: qsa_init_state() unimplemented!\n");
}

void evolution_szegedy(const double& t, const int& n){
    (void)t,(void)n;
    //TODO: implement
    throw std::runtime_error("ERROR: evolution_szegedy() unimplemented!\n");
}

void evolution_measure(const double& t, const int& n){
    (void)t,(void)n;
    //TODO: implement
    throw std::runtime_error("ERROR: evolution_measure() unimplemented!\n");
}

void evolution_tracing(const double& t, const int& n){
    (void)t,(void)n;
    //TODO: implement
    throw std::runtime_error("ERROR: evolution_tracing() unimplemented!\n");
//    (void)n;
//  for (uint i = 0; i < 3; i++) {
//    suqa::apply_pauli_TP_rotation({bm_spin_tilde[(0+i)%3],bm_spin_tilde[(1+i)%3]}, {PAULI_X,PAULI_X}, -t);
//  }
//
}

/* Measure facilities */
const uint op_bits = 3; // 2^op_bits is the number of eigenvalues for the observable
const bmReg bm_op = bm_q; // where the measure has to be taken
const std::vector<double> op_vals = {2.0,0.0,-2.0, 0.0,0.0,0.0,0.0,0.0}; // eigvals

 
// change basis to the observable basis somewhere in the system registers
void apply_measure_rotation(){
//    self_plaquette(bm_qlink1, bm_qlink0, bm_qlink2, bm_qlink0);
}

// inverse of the above function
void apply_measure_antirotation(){
//    inverse_self_plaquette(bm_qlink1, bm_qlink0, bm_qlink2, bm_qlink0);
}

// map the classical measure recorded in creg_vals
// to the corresponding value of the observable;
// there is no need to change it
double get_meas_opvals(const uint& creg_vals){
    return op_vals[creg_vals];
}

// actually perform the measure
// there is no need to change it
double measure_X(pcg& rgen){
    std::vector<uint> classics(op_bits);
    
    apply_measure_rotation();

    std::vector<double> rdoubs(op_bits);
    for(auto& el : rdoubs){
        el = rgen.doub();
    }
    suqa::measure_qbits(bm_op, classics, rdoubs);

    apply_measure_antirotation();

    uint meas = 0U;
    for(uint i=0; i<op_bits; ++i){
        meas |= (classics[i] << i);
    }

    return get_meas_opvals(meas);
}

/* Moves facilities */
#define NMoves 5

//std::vector<double> C_weightsums(NMoves);

std::vector<double> C_weightsums;


void apply_C(const uint &Ci,double rot_angle){
    // move 0 -> Ci=0, inverse move 0 -> Ci=9
    //bool is_inverse = Ci>=HNMoves;
    //double actual_angle = (is_inverse)? -rot_angle : rot_angle;
//    double actual_angle = rot_angle;
    switch (Ci){
        case 0:
        case 1:
        case 2:
            suqa::apply_h(Ci);
            break;
        case 3:
        case 4:
            suqa::apply_cx(Ci-3,Ci-2);
            break;

        default:
            throw std::runtime_error("ERROR: apply_C() unimplemented!\n");
    }


}

void apply_C_inverse(const uint &Ci,double rot_angle){
    apply_C(Ci,-rot_angle);
    // or, equivalent:
//    apply_C((Ci+HNMoves)%NMoves,rot_angle);
//    throw std::runtime_error("ERROR: apply_C_inverse() unimplemented!\n");
}

void qsa_apply_C(const uint &Ci){
    (void)Ci;
    //TODO: implement
    throw std::runtime_error("ERROR: qsa_apply_C() unimplemented!\n");
//  suqa::apply_h(bm_spin_tilde[Ci]);
// suqa::apply_h(state,bm_spin_tilde[(Ci+1)%3]);


  // suqa::apply_h(state,bm_spin_tilde);
}

void qsa_apply_C_inverse(const uint &Ci){
    (void)Ci;
    //TODO: implement
    throw std::runtime_error("ERROR: qsa_apply_C() unimplemented!\n");
//  if(Ci>2) throw std::runtime_error("ERROR: wrong move selection");
//  //suqa::apply_h(state,bm_spin_tilde);
//  //suqa::apply_h(state,bm_spin_tilde[(Ci+1)%3]);
//  suqa::apply_h(bm_spin_tilde[Ci]);
}

//std::vector<double> get_C_weightsums(){ return C_weightsums; }

std::vector<double> get_C_weightsums(){ 
    static bool init_done=false;
    // first initialization
    if(not init_done){
        C_weightsums=std::vector<double>(NMoves);
        for(int i=1; i<=NMoves; ++i){
            C_weightsums[i-1]=i/(double)NMoves;
        }
        init_done=true;
    }
    return C_weightsums; 
}

