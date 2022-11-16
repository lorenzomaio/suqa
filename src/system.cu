#include "system.cuh"

/* d4 gauge theory - one plaquette
 
    link state 3 qubits
    system state: 2 links -> 6 qubits
    +1 ancillary qubit

    .   .  
    1   
    o 0 .  

    operation table for the D4 group:

    [abc -> a*4+b*2+c]

    representation: ρ_{abc} = (-1)^b [[0,1],[1,0]]^a [[i, 0],[0,-i]]^c

    ρ_{a'b'c'}ρ_{abc} = (-1)^{b'+b} 

     
          0   1   2   3   4   5   6   7
        ________________________________
       |
    0  |  0   1   2   3   4   5   6   7
       |
    1  |  1   2   3   0   7   4   5   6 
       | 
    2  |  2   3   0   1   6   7   4   5
       |
    3  |  3   0   1   2   5   6   7   4
       |
    4  |  4   7   6   5   0   1   2   3
       |
    5  |  5   4   7   6   1   0   1   2
       |
    6  |  6   5   4   7   2   1   0   1
       |
    7  |  7   6   5   4   3   2   1   0

 */



double g_beta;

__inline__ double f1(double b){
    return log((3+cosh(2.*b))/(2*sinh(b)*sinh(b)));
}

__inline__ double f2(double b){
    return -log(tanh(b));
}

void init_state(){
    //suqa::init_state();

    //suqa::apply_h(bm_qlink0[0]);
    //suqa::apply_cx(bm_qlink0[0], bm_qlink1[0]);
    //suqa::apply_h(bm_qlink0[1]);
    //suqa::apply_cx(bm_qlink0[1], bm_qlink1[1]);
    //suqa::apply_h(bm_qlink0[2]);
    //suqa::apply_cx(bm_qlink0[2], bm_qlink1[2]);
    //suqa::apply_mcx({bm_qlink1[0], bm_qlink1[2]}, {0U,1U}, bm_qlink1[1]);

    suqa::init_state();

    suqa::apply_h(bm_qlink0[0]);
    suqa::apply_h(bm_qlink0[1]);
    suqa::apply_h(bm_qlink0[2]);
    suqa::apply_h(bm_qlink1[0]);
    suqa::apply_h(bm_qlink1[1]);
    suqa::apply_h(bm_qlink1[2]);
    
// automatically gauge-invariant initialization
/*    .   .   .
 *    1   2
 *    o 0 o I .
 *   g1  g2
 *  
 *  U0 -> g2 U0 g1' -> g
 *  U1 -> g1 U1 g1' -> I
 *  U2 -> g2 U2 g2' -> I
 *  U3 -> g1 U3 g2' -> g'
 *
 *  g=g2 g1', sum over g
 *
 *  GF: g2=g1 U3
 *  
 *  U0 -> U0n = g1 U3 U0 g1'     -> I
 *  U1 -> U1n = g1 U1 g1'        -> I
 *  U2 -> U2n = g1 U3 U2 U3' g1' -> I
 *  U3 -> U3n = I                -> I
 *
 */

}


/* Quantum evolutor of the state */

void inversion(const bmReg& q){
    suqa::apply_mcx({q[0],q[2]},{1U,0U},q[1]); 
}

void left_multiplication(const bmReg& qr1, const bmReg& qr2){
    // applies group element from register qr1 to register qr2
    // |...,U_{qr1},...,U_{qr2},...> -> |...,U_{qr1},...,U_{qr1}U_{qr2},...>
    suqa::apply_cx(qr1[1], qr2[1]);
    suqa::apply_mcx({qr1[0], qr2[0]}, qr2[1]);
    suqa::apply_cx(qr1[0], qr2[0]);
    suqa::apply_mcx({qr1[0], qr2[2]}, qr2[1]);
    suqa::apply_cx(qr1[2], qr2[2]);
}

void self_plaquette(const bmReg& qr0, const bmReg& qr1, const bmReg& qr2, const bmReg& qr3){
    // applies the following operation in the computational basis
    // |...,U_{qr0}...,U_{qr1}...,U_{qr2},...,U_{qr3},...> -> 
    //                      .
    //                      .
    //                      V
    // |...,U_{qr3}U'_{qr2}U'_{qr1}U_{qr0}...,U_{qr1}...,U_{qr2},...,U_{qr3},...>
    inversion(qr1);
    left_multiplication(qr1, qr0);
    inversion(qr1);
    inversion(qr2);
    left_multiplication(qr2, qr0);
    inversion(qr2);
    left_multiplication(qr3, qr0);
}

void inverse_self_plaquette(const bmReg& qr0, const bmReg& qr1, const bmReg& qr2, const bmReg& qr3){
    // inverse operation of self_plaquette
    inversion(qr3);
    left_multiplication(qr3, qr0);
    inversion(qr3);
    left_multiplication(qr2, qr0);
    left_multiplication(qr1, qr0);
}



void cphases(uint qaux, uint q0b, double alpha1, double alpha2){
    // eigenvalues of the trace operator
    suqa::apply_cx(qaux, q0b);
    suqa::apply_cu1(q0b, qaux, alpha1, 1U);
    suqa::apply_cx(qaux, q0b);
    suqa::apply_cu1(q0b, qaux, alpha2, 1U);
}

void self_trace_operator(const bmReg& qr, const bmReg& qr1, double th){
    //suqa::apply_mcx({qr[0],qr[2]}, {0U,0U}, qaux); 
    //cphases(qaux, qr[1], th, -th);
    //suqa::apply_mcx({qr[0],qr[2]}, {0U,0U}, qaux); 

    suqa::apply_mcu1({qr[2],qr1[2]}, {0U,0U}, qr[1], th);
    suqa::apply_mcu1({qr[2],qr1[2], qr[0]}, {0U,1U,0U}, qr[1], th);
    suqa::apply_mcu1({qr[2],qr1[2], qr[0]}, {0U,1U,1U}, qr[1], -th);
    suqa::apply_mcu1({qr[2],qr1[2], qr1[0]}, {1U,0U,0U}, qr[1], th);
    suqa::apply_mcu1({qr[2],qr1[2], qr1[0]}, {1U,0U,1U}, qr[1], -th);
    suqa::apply_mcu1({qr[2],qr1[2], qr[0],qr1[0]}, {1U,1U,0U,0U}, qr[1], th);
    suqa::apply_mcu1({qr[2],qr1[2], qr[0],qr1[0]}, {1U,1U,1U,1U}, qr[1], th);
    suqa::apply_mcu1({qr[2],qr1[2], qr[0],qr1[0]}, {1U,1U,0U,1U}, qr[1], -th);
    suqa::apply_mcu1({qr[2],qr1[2], qr[0],qr1[0]}, {1U,1U,1U,0U}, qr[1], -th);   
    
// Alternative implementation (of the previous implementation)
//    suqa::apply_mcu1({qr[0],qr[2]}, {0U,0U}, qr[1],th);  // u1(θ) = [[1,0],[0,e^{iθ}]]
//    suqa::apply_mcx({qr[0],qr[2]}, {0U,0U}, qr[1]);  
//    suqa::apply_mcu1({qr[0],qr[2]}, {0U,0U}, qr[1],-th);
//    suqa::apply_mcx({qr[0],qr[2]}, {0U,0U}, qr[1]);
}

void fourier_transf_d4(const bmReg& qr){
    suqa::apply_cx(qr[2], qr[0]);
    suqa::apply_cx(qr[0], qr[2]);
    suqa::apply_tdg(qr[2]);
    suqa::apply_tdg(qr[2]);
    suqa::apply_cx(qr[1], qr[2]);
    suqa::apply_h(qr[0]);
    suqa::apply_h(qr[1]);
    suqa::apply_h(qr[2]);
    suqa::apply_t(qr[1]);
    suqa::apply_tdg(qr[2]);
    suqa::apply_cx(qr[1], qr[2]);
    suqa::apply_cx(qr[0], qr[1]);
    suqa::apply_h(qr[1]);
    suqa::apply_t(qr[1]);
    suqa::apply_t(qr[1]);
    suqa::apply_h(qr[1]);
}


void inverse_fourier_transf_d4(const bmReg& qr){
    suqa::apply_h(qr[1]);
    suqa::apply_tdg(qr[1]);
    suqa::apply_tdg(qr[1]);
    suqa::apply_h(qr[1]);
    suqa::apply_cx(qr[0], qr[1]);
    suqa::apply_cx(qr[1], qr[2]);
    suqa::apply_t(qr[2]);
    suqa::apply_tdg(qr[1]);
    suqa::apply_h(qr[0]);
    suqa::apply_h(qr[1]);
    suqa::apply_h(qr[2]);
    suqa::apply_cx(qr[1], qr[2]);
    suqa::apply_t(qr[2]);
    suqa::apply_t(qr[2]);
    suqa::apply_cx(qr[0], qr[2]);
    suqa::apply_cx(qr[2], qr[0]);
}

void momentum_phase(const bmReg& qr, const uint& qaux, double th1, double th2){
    suqa::apply_mcx(qr, {0U,0U,0U}, qaux);
    DEBUG_CALL(printf("\tafter suqa::apply_mcx(qr, {0U,0U,0U}, qaux)\n"));
    DEBUG_READ_STATE();
    suqa::apply_cx(qaux, qr[2]);
    DEBUG_CALL(printf("\tafter suqa::apply_cx(qaux=%u, qr[2]=%u)\n",qaux,qr[2]));
    DEBUG_READ_STATE();
    suqa::apply_cu1(qaux, qr[2], th1);
    DEBUG_CALL(printf("\tafter suqa::apply_cu1(qaux, qr[2],%lg)\n",th1));
    DEBUG_READ_STATE();
    suqa::apply_cx(qaux, qr[2]);
    DEBUG_CALL(printf("\tafter suqa::apply_cu1(qaux, qr[2], th1)\n"));
    DEBUG_READ_STATE();
    suqa::apply_mcx(qr, {0U,0U,0U}, qaux);
    DEBUG_CALL(printf("\tafter suqa::apply_mcx(qr, {0U,0U,0U}, qaux)\n"));
    DEBUG_READ_STATE();


    // Alternative implementation without qaux
//    suqa::apply_mcx(qr, {0U,0U,0U}, qr[2]);
//    suqa::apply_mcu1(qr, {0U,0U,0U}, qr[2],th1);
//    suqa::apply_mcx(qr, {0U,0U,0U}, qr[2]);

    suqa::apply_u1(qr[2], th2);
    DEBUG_CALL(printf("\tafter suqa::apply_u1(qr[2], th2)\n"));
    DEBUG_READ_STATE();
}

void evolution(const double& t, const int& n){
    const double dt = t/(double)n;

    const double theta1 = dt*f1(g_beta);    // eigenvalues of kinetic hamiltonian on single gauge variable
    const double theta2 = dt*f2(g_beta);
    const double theta = 2*dt*g_beta;       // see Lamm's paper (the factor 2 is included here)
    DEBUG_CALL(printf("Evolution parameters:\ng_beta = %.16lg, dt = %.16lg, thetas: %.16lg %.16lg %.16lg\n", g_beta, dt, theta1, theta2, theta));

    for(uint ti=0; ti<(uint)n; ++ti){
        // V^{1/2}

        self_trace_operator(bm_qlink1, bm_qlink0, theta);
        DEBUG_CALL(printf("after self_trace_operator()\n"));
        DEBUG_READ_STATE();

        // T
        fourier_transf_d4(bm_qlink0);
        DEBUG_CALL(printf("after fourier_transf_d4(bm_qlink0)\n"));
        DEBUG_READ_STATE();
        fourier_transf_d4(bm_qlink1);
        DEBUG_CALL(printf("after fourier_transf_d4(bm_qlink1)\n"));
        DEBUG_READ_STATE();

        momentum_phase(bm_qlink0, bm_qaux[0], theta1, theta2);
        DEBUG_CALL(printf("after momentum_phase(bm_qlink0, bm_qaux[0], theta1, theta2)\n"));
        DEBUG_READ_STATE();
        momentum_phase(bm_qlink1, bm_qaux[0], theta1, theta2);
        DEBUG_CALL(printf("after momentum_phase(bm_qlink1, bm_qaux[0], theta1, theta2)\n"));
        DEBUG_READ_STATE();

        inverse_fourier_transf_d4(bm_qlink1);
        DEBUG_CALL(printf("after inverse_fourier_transf_d4(bm_qlink1)\n"));
        DEBUG_READ_STATE();
        inverse_fourier_transf_d4(bm_qlink0);
        DEBUG_CALL(printf("after inverse_fourier_transf_d4(bm_qlink0)\n"));
        DEBUG_READ_STATE();
//
//
//        // V^{1/2}
//
//        self_plaquette(bm_qlink2, bm_qlink3, bm_qlink1, bm_qlink3);
//        DEBUG_CALL(printf("after self_plaquette2()\n"));
//        DEBUG_READ_STATE();
//        self_trace_operator(bm_qlink2, bm_qaux[0], theta*0.5);
//        DEBUG_CALL(printf("after self_trace_operator()\n"));
//        DEBUG_READ_STATE();
//        inverse_self_plaquette(bm_qlink2, bm_qlink3, bm_qlink1, bm_qlink3);
//        DEBUG_CALL(printf("after inverse_self_plaquette2()\n"));
//        DEBUG_READ_STATE();
//
//        self_plaquette(bm_qlink1, bm_qlink0, bm_qlink2, bm_qlink0);
//        DEBUG_CALL(printf("after self_plaquette1()\n"));
//        DEBUG_READ_STATE();
//        self_trace_operator(bm_qlink1, bm_qaux[0], theta*0.5);
//        DEBUG_CALL(printf("after self_trace_operator()\n"));
//        DEBUG_READ_STATE();
//        inverse_self_plaquette(bm_qlink1, bm_qlink0, bm_qlink2, bm_qlink0);
//        DEBUG_CALL(printf("after inverse_self_plaquette1()\n"));
//        DEBUG_READ_STATE();
    }
}



// qsa specifics
void qsa_init_state(){
    //TODO: implement
    throw std::runtime_error("ERROR: qsa_init_state() unimplemented!\n");
//    suqa::init_state();
//    suqa::apply_h(bm_spin[0]);
//    suqa::apply_h(bm_spin[1]);
//    suqa::apply_h(bm_spin[2]);
//    suqa::apply_cx(bm_spin[0], bm_spin_tilde[0]);
//    suqa::apply_cx(bm_spin[1], bm_spin_tilde[1]);
//    suqa::apply_cx(bm_spin[2], bm_spin_tilde[2]);
}

void evolution_szegedy(const double& t, const int& n){
    (void)t,(void)n;
    //TODO: implement
    throw std::runtime_error("ERROR: evolution_szegedy() unimplemented!\n");
//    (void)n;
//      DEBUG_CALL(std::cout<<"before evolution_szegedy()"<<std::endl);
//      DEBUG_READ_STATE();
//      DEBUG_CALL(std::cout<<"apply evolution_szegedy()"<<std::endl);
//    for (uint i = 0; i < 3; i++) {
//      suqa::apply_pauli_TP_rotation({bm_spin_tilde[(0+i)%3],bm_spin_tilde[(1+i)%3]}, {PAULI_X,PAULI_X}, -t);
//      DEBUG_CALL(std::cout<<"apply pauli 1 it "<<i<<std::endl);
//      DEBUG_READ_STATE();
//      suqa::apply_pauli_TP_rotation({bm_spin[(0+i)%3],bm_spin[(1+i)%3]}, {PAULI_X,PAULI_X}, t);
//      DEBUG_CALL(std::cout<<"apply pauli 2 it "<<i<<std::endl);
//      DEBUG_READ_STATE();
//    }
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
const bmReg bm_op = bm_qlink1; // where the measure has to be taken
const std::vector<double> op_vals = {2.0,0.0,-2.0, 0.0,0.0,0.0,0.0,0.0}; // eigvals

 
// change basis to the observable basis somewhere in the system registers
void apply_measure_rotation(){
    //self_plaquette(bm_qlink1, bm_qlink0, bm_qlink2, bm_qlink0);
    //self_trace_operator(bm_qlink1, bm_qlink0, theta);
}

// inverse of the above function
void apply_measure_antirotation(){
    //self_trace_operator(bm_qlink1, bm_qlink0, -theta);
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
    
    //apply_measure_rotation();

    std::vector<double> rdoubs(op_bits);
    for(auto& el : rdoubs){
        el = rgen.doub();
    }
    suqa::measure_qbits(bm_op, classics, rdoubs);

    //apply_measure_antirotation();

    uint meas = 0U;
    for(uint i=0; i<op_bits; ++i){
        meas |= (classics[i] << i);
    }

    return get_meas_opvals(meas);
}

/* Moves facilities */
#define NMoves 3

//std::vector<double> C_weightsums(NMoves);

std::vector<double> C_weightsums = {1./3., 2./3., 1.0};
//#define HNMoves (NMoves>>1)


void apply_C(const uint &Ci,double rot_angle){
    // move 0 -> Ci=0, inverse move 0 -> Ci=9
    //bool is_inverse = Ci>=HNMoves;
    //double actual_angle = (is_inverse)? -rot_angle : rot_angle;
    double actual_angle = rot_angle;
    switch (Ci){
        case 0:
        case 1:// eigenvalues of kinetic hamiltonian on single gauge variable
        {
            const double theta1 = actual_angle*f1(g_beta);    
            const double theta2 = actual_angle*f2(g_beta);
            DEBUG_CALL(printf("actual_angle = %lg; g_beta = %lg; f1 = %lg; f2 = %lg\n",actual_angle,g_beta,f1(g_beta), f2(g_beta)));
            fourier_transf_d4(bm_qlinks[Ci]);
            momentum_phase(bm_qlinks[Ci], bm_qaux[0], theta1, theta2);
            inverse_fourier_transf_d4(bm_qlinks[Ci]);
            break;
        }
        case 2: // plaquette
        {
            self_trace_operator(bm_qlink1, bm_qlink0, actual_angle);
            break;
        }
        
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

std::vector<double> get_C_weightsums(){ return C_weightsums; }

// std::vector<double> get_C_weightsums(){ 
//    static bool init_done=false;
    // first initialization
//    if(not init_done){
//        for(int i=1; i<=NMoves; ++i){
//            C_weightsums[i-1]=i/(double)NMoves;
//        }
//        init_done=true;
//    }
//    return C_weightsums; }

