#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include "sine.h"
#include "cosine.h"
#include "inv_matrix.h"

// reference: [1]Liu JinKun. Robot Control System Design and MATLAB Simulation[M]. Tsinghua University Press, 2008.
// [2]Shuzhi S G, Hang C C, Woon L C. Adaptive neural network control of robot manipulators in task space[J]. IEEE transactions on industrial electronics, 1997, 44(6): 746-752.

// global variables declaration
#define PI 3.14159
#define H 7              // input layer neurons number
#define IN 4             // hidden layer neurons number
#define OUT 2            // output layer neurons number
#define ARRAY_SIZE 10000 // sampling times

static double Ts = 0.001;                                                                                                                                                                 // sampling period
static double t0 = 0.0;                                                                                                                                                                   // start time
static double t1 = 10.0;                                                                                                                                                                  // end time
static double center_D[OUT][H] = {{0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75}, {0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75}};                                                                          // inertia matrix RBF function center
static double center_G[OUT][H] = {{0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75}, {0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75}};                                                                          // gravity matrix RBF function center
static double center_C[IN][H] = {{0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75}, {0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75}, {-1.5, -1, -0.5, 0, 0.5, 1.0, 1.50}, {-1.5, -1, -0.5, 0, 0.5, 1.0, 1.50}}; // Coriolis force and centrifugal force matrix RBF function center
static double width = 10;                                                                                                                                                                 // RBF function width

double phi_D11[H], phi_D12[H], phi_D21[H], phi_D22[H];
double phi_G1[H], phi_G2[H];
double phi_C11[H], phi_C12[H], phi_C21[H], phi_C22[H];
double weight_D11[H], weight_D12[H], weight_D21[H], weight_D22[H];
double weight_G1[H], weight_G2[H];
double weight_C11[H], weight_C12[H], weight_C21[H], weight_C22[H];
double derivative_weight_D11[H], derivative_weight_D12[H], derivative_weight_D21[H], derivative_weight_D22[H];
double derivative_weight_G1[H], derivative_weight_G2[H];
double derivative_weight_C11[H], derivative_weight_C12[H], derivative_weight_C21[H], derivative_weight_C22[H];

void matrix_Multi(double C[][2], double A[][2], double B[][2], int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        for (int k = 0; k < cols2; k++){
            C[j][k] = 0.0;
            for (int g = 0; g < cols1; g++){
                C[j][k] += A[j][g] * B[g][k];
            }
        }
    }
}

void matrix_Multi1(double C[], double A[][2], double B[], int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        C[j] = 0.0;
        for (int k = 0; k < cols2; k++){
            C[j] += A[j][k] * B[k];
        }
    }
}

struct _archive{
    double x1_archive[ARRAY_SIZE];
    double dx1_archive[ARRAY_SIZE];
    double ddx1_archive[ARRAY_SIZE];
    double x2_archive[ARRAY_SIZE];
    double dx2_archive[ARRAY_SIZE];
    double ddx2_archive[ARRAY_SIZE];
    double error1_archive[ARRAY_SIZE];
    double error2_archive[ARRAY_SIZE];
    double error1_velocity_archive[ARRAY_SIZE];
    double error2_velocity_archive[ARRAY_SIZE];
    double Fx1_archive[ARRAY_SIZE];
    double Fx2_archive[ARRAY_SIZE];
    double tol1_archive[ARRAY_SIZE];
    double tol2_archive[ARRAY_SIZE];
    double D_task_estimate_norm_archive[ARRAY_SIZE];
    double D_task_norm_archive[ARRAY_SIZE];
    double G_task_estimate_norm_archive[ARRAY_SIZE];
    double G_task_norm_archive[ARRAY_SIZE];
    double C_task_estimate_norm_archive[ARRAY_SIZE];
    double C_task_norm_archive[ARRAY_SIZE];
} archive;

Data xd1, dxd1, ddxd1, xd2, dxd2, ddxd2;

struct Amp{
    double xd1;
    double dxd1;
    double ddxd1;
    double xd2;
    double dxd2;
    double ddxd2;
};

struct M0{
    double xd1;
    double dxd1;
    double ddxd1;
    double xd2;
    double dxd2;
    double ddxd2;
};

struct B0{
    double xd1;
    double dxd1;
    double ddxd1;
    double xd2;
    double dxd2;
    double ddxd2;
};

void Input(Data *xd1, Data *dxd1, Data *ddxd1, Data *xd2, Data *dxd2, Data *ddxd2, double Ts, double t0, double t1){

    struct Amp amp; // amplitude
    amp.xd1 = 0.2;
    amp.dxd1 = -0.2 * PI;
    amp.ddxd1 = -0.2 * pow(PI, 2);
    amp.xd2 = 0.2;
    amp.dxd2 = 0.2 * PI;
    amp.ddxd2 = -0.2 * pow(PI, 2);

    struct M0 m0; // angular frequency
    m0.xd1 = PI;
    m0.dxd1 = PI;
    m0.ddxd1 = PI;
    m0.xd2 = PI;
    m0.dxd2 = PI;
    m0.ddxd2 = PI;

    struct B0 b0; // vertical shift
    b0.xd1 = 1.0;
    b0.dxd1 = 0.0;
    b0.ddxd1 = 0.0;
    b0.xd2 = 1.0;
    b0.dxd2 = 0.0;
    b0.ddxd2 = 0.0;

    cosine(xd1, Ts, t0, t1, amp.xd1, m0.xd1, b0.xd1);         // desired Cartesian coordinate displacement of link 1
    sine(dxd1, Ts, t0, t1, amp.dxd1, m0.dxd1, b0.dxd1);       // desired Cartesian coordinate velocity of link 1
    cosine(ddxd1, Ts, t0, t1, amp.ddxd1, m0.ddxd1, b0.ddxd1); // desired Cartesian coordinate acceleration of link 1
    sine(xd2, Ts, t0, t1, amp.xd2, m0.xd2, b0.xd2);           // desired Cartesian coordinate displacement of link 2
    cosine(dxd2, Ts, t0, t1, amp.dxd2, m0.dxd2, b0.dxd2);     // desired Cartesian coordinate velocity of link 2
    sine(ddxd2, Ts, t0, t1, amp.ddxd2, m0.ddxd2, b0.ddxd2);   // desired angular acceleration of link 2
}

struct _system_state{
    double x1;   // actual Cartesian coordinate displacement of link 1
    double dx1;  // actual Cartesian coordinate velocity of link 1
    double ddx1; // actual Cartesian coordinate acceleration of link 1
    double x2;   // actual Cartesian coordinate displacement of link 2
    double dx2;  // actual Cartesian coordinate velocity of link 2
    double ddx2; // actual Cartesian coordinate acceleration of link 2
    double q1;   // actual angular displacement of link 1
    double dq1;  // actual angular velocity of link 1
    double q2;   // actual angular displacement of link 2
    double dq2;  // actual angular velocity of link 2
    double p1;   // actual absolute angular displacement of link 1
    double dp1;  // actual absolute angular velocity of link 1
    double p2;   // actual absolute angular displacement of link 2
    double dp2;  // actual absolute angular velocity of link 2
} system_state;

struct _control_variable{
    double force1; // control force of link 1
    double force2; // control force of link 2
    double tol1;   // control torque of link 1
    double tol2;   // control torque of link 2
} control_variable;

struct _parameter{
    double l1;   // control force of link 1
    double l2;   // control force of link 2
    double g;    // control torque of link 1
    double P[H]; // M = P + pl * L
    double pl;
    double M[H]; // matrix parameter mi in joint control mi
    double L[H];
} parameter;

struct _dynamics{
    double D_joint[OUT][OUT];             // inertia matrix for model joint space
    double G_joint[OUT];                  // gravity matrix for model joint space
    double C_joint[OUT][OUT];             // Coriolis matrix for model joint space
    double D_task[OUT][OUT];              // inertia matrix for model task space
    double G_task[OUT];                   // gravity matrix for model task space
    double C_task[OUT][OUT];              // Coriolis matrix for model task space
    double D_task_norm;                   // two-paradigm number of D_task
    double G_task_norm;                   // two-paradigm number of G_task
    double C_task_norm;                   // two-paradigm number of C_task
    double DSNN_estimate[OUT][OUT];       // estimate of RBF network modeling term DSNN
    double GSNN_estimate[OUT];            // estimate of RBF network modeling term GSNN
    double CDNN_estimate[OUT][OUT];       // estimate of RBF network modeling term CDNN
    double D_estimate_norm;               // second-paradigm number of estimate of RBF network modeling term DSNN
    double G_estimate_norm;               // second-paradigm number of estimate of RBF network modeling term GSNN
    double C_estimate_norm;               // second-paradigm number of estimate of RBF network modeling term CDNN
    double inv_D_task[OUT][OUT];          // inverse of inertia matrix for model task space
    double Jacobian[OUT][OUT];            // Jacobi matrix
    double Jacobian_derivative[OUT][OUT]; // derivative of Jacobi matrix
} dynamics;

struct _controller{
    double controller_u1;
    double controller_u2;
    double controller_u3;
    double controller_u4;
    double controller_u5;
    double controller_u6;
    double controller_u7;
    double controller_u8;
    double controller_u9;
    double controller_u10;
    double controller_out1;
    double controller_out2;
    double controller_out3;
    double controller_out4;
    double controller_out5;
    double err1;                                                       // Cartesian coordinate displacement error of link 1
    double err1_velocity;                                              // Cartesian coordinate velocity error of link 1
    double err2;                                                       // Cartesian coordinate displacement error of link 2
    double err2_velocity;                                              // Cartesian coordinate velocity error of link 2
    double Tau_D11[H][H], Tau_D12[H][H], Tau_D21[H][H], Tau_D22[H][H]; // Eq. 3.60 and 3.77 define, adaptive law design
    double Tau_G1[H][H], Tau_G2[H][H];                                 // Eq. 3.62 and 3.77 define, adaptive law design
    double Tau_C11[H][H], Tau_C12[H][H], Tau_C21[H][H], Tau_C22[H][H]; // Eq. 3.61 and 3.77 define, adaptive law design
    double Lambda[OUT][OUT];                                           // error's weight factor
    double r1;                                                         // refer to 3.5.3 controller design
    double r2;
    double dxr1;  // derivative of xr1
    double dxr2;  // derivative of xr2
    double ddxr1; // second-order derivative of xr1
    double ddxr2; // second-order derivative of xr2
    double K1;    // controller gain of link 1
    double K2;    // controller gain of link 2
    double Ks;    // robust term gain
} controller;

void CONTROL_init(){
    system_state.x1 = 1.0;
    system_state.dx1 = 0.0;
    system_state.x2 = 1.0;
    system_state.dx2 = 0.0;
    controller.controller_u1 = xd1.y[0];
    controller.controller_u2 = dxd1.y[0];
    controller.controller_u3 = ddxd1.y[0];
    controller.controller_u4 = xd2.y[0];
    controller.controller_u5 = dxd2.y[0];
    controller.controller_u6 = ddxd2.y[0];
    controller.controller_u7 = system_state.x1;
    controller.controller_u8 = system_state.dx1;
    controller.controller_u9 = system_state.x2;
    controller.controller_u10 = system_state.dx2;
    controller.err1 = xd1.y[0] - system_state.x1;
    controller.err1_velocity = dxd1.y[0] - system_state.dx1;
    controller.err2 = xd1.y[0] - system_state.x2;
    controller.err2_velocity = dxd2.y[0] - system_state.dx2;

    for (int j = 0; j < H; j++){
        for (int k = 0; k < H; k++){
            if (j == k){
                controller.Tau_D11[j][k] = 2.0;
                controller.Tau_D12[j][k] = 2.0;
                controller.Tau_D21[j][k] = 2.0;
                controller.Tau_D22[j][k] = 2.0;
            }
            else{
                controller.Tau_D11[j][k] = 0.0;
                controller.Tau_D12[j][k] = 0.0;
                controller.Tau_D21[j][k] = 0.0;
                controller.Tau_D22[j][k] = 0.0;
            }
        }
    }

    for (int j = 0; j < H; j++){
        for (int k = 0; k < H; k++){
            if (j == k){
                controller.Tau_G1[j][k] = 5.0;
                controller.Tau_G2[j][k] = 5.0;
            }
            else{
                controller.Tau_G1[j][k] = 0.0;
                controller.Tau_G2[j][k] = 0.0;
            }
        }
    }

    for (int j = 0; j < H; j++){
        for (int k = 0; k < H; k++){
            if (j == k){
                controller.Tau_C11[j][k] = 0.5;
                controller.Tau_C12[j][k] = 0.5;
                controller.Tau_C21[j][k] = 0.5;
                controller.Tau_C22[j][k] = 0.5;
            }
            else{
                controller.Tau_C11[j][k] = 0.0;
                controller.Tau_C12[j][k] = 0.0;
                controller.Tau_C21[j][k] = 0.0;
                controller.Tau_C22[j][k] = 0.0;
            }
        }
    }

    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            if (j == k){
                controller.Lambda[j][k] = 15.0;
            }
            else{
                controller.Lambda[j][k] = 0.0;
            }
        }
    }

    controller.r1 = controller.err1_velocity + controller.Lambda[0][0] * controller.err1;
    controller.r2 = controller.err2_velocity + controller.Lambda[1][1] * controller.err2;
    controller.K1 = 50;
    controller.K2 = 40;
    controller.Ks = 0.6;
}

struct _plant{
    double plant_u1;
    double plant_u2;
    double plant_out1;
    double plant_out2;
    double plant_out3;
    double plant_out4;
    double plant_out5;
    double plant_out6;
    double plant_out7;
    double plant_out8;
    double plant_out9;
} plant;

void PLANT_init(){
    system_state.x1 = 1.0;
    system_state.dx1 = 0.0;
    system_state.x2 = 1.0;
    system_state.dx2 = 0.0;
    plant.plant_u1 = 0.0;
    plant.plant_u2 = 0.0;
    plant.plant_out1 = system_state.x1;
    plant.plant_out2 = system_state.dx1;
    plant.plant_out3 = system_state.x2;
    plant.plant_out4 = system_state.dx2;
}

double PLANT_realize(int i){
    plant.plant_u1 = control_variable.force1;
    plant.plant_u2 = control_variable.force2;

    double P_value[] = {1.66, 0.42, 0.63, 3.75, 1.25};
    for (int j = 0; j < H; ++j){
        parameter.P[j] = P_value[j];
    }
    parameter.g = 9.8;
    parameter.l1 = 1;
    parameter.l2 = 1;
    double L_value[] = {pow(parameter.l1, 2), pow(parameter.l2, 2), parameter.l1 * parameter.l2, parameter.l1, parameter.l2};
    for (int j = 0; j < H; ++j){
        parameter.L[j] = L_value[j];
    }

    if (i * Ts > 4.0){
        parameter.pl = 0.0;
    }
    else{
        parameter.pl = 0.5;
    }

    for (int j = 0; j < H; ++j){
        parameter.M[j] = parameter.P[j] + parameter.pl * parameter.L[j];
    }

    double Q = (pow(system_state.x1, 2) + pow(system_state.x2, 2) - pow(parameter.l1, 2) - pow(parameter.l2, 2)) / (2 * parameter.l1 * sqrt(pow(system_state.x1, 2) + pow(system_state.x2, 2)));
    system_state.q2 = acos(Q); // Eq. 3.67
    system_state.dq2 = -1 / sqrt(1 - pow(Q, 2));
    double ratio = system_state.x2 / system_state.x1;
    system_state.p1 = atan(ratio);
    system_state.dp1 = 1 / (1 + pow(ratio, 2));
    double B = sqrt(pow(system_state.x1, 2) + pow(system_state.x2, 2) + pow(parameter.l1, 2) - pow(parameter.l2, 2)) / (2 * parameter.l1 * sqrt(pow(system_state.x1, 2) + pow(system_state.x2, 2)));
    system_state.p2 = acos(B);
    system_state.dp2 = -1 / sqrt(1 - pow(B, 2));

    if (system_state.q2 > 0){
        system_state.q1 = system_state.p1 - system_state.p2;
        system_state.dq1 = system_state.dp1 - system_state.dp2;
    }
    else{
        system_state.q1 = system_state.p1 + system_state.p2;
        system_state.dq1 = system_state.dp1 + system_state.dp2;
    }

    dynamics.Jacobian[0][0] = -parameter.l1 * (system_state.q1) - parameter.l2 * sin(system_state.q1 + system_state.q2); // Jacobi matrix
    dynamics.Jacobian[0][1] = -parameter.l2 * sin(system_state.q1 + system_state.q2);
    dynamics.Jacobian[1][0] = parameter.l1 * cos(system_state.q1) + parameter.l2 * cos(system_state.q1 + system_state.q2);
    dynamics.Jacobian[1][1] = parameter.l2 * cos(system_state.q1 + system_state.q2);

    dynamics.Jacobian_derivative[0][0] = -system_state.dq1 * cos(system_state.q1) - (system_state.dq1 + system_state.dq2) * cos(system_state.q1 + system_state.q2); // derivative of Jacobi matrix
    dynamics.Jacobian_derivative[0][1] = -(system_state.dq1 + system_state.dq2) * cos(system_state.q1 + system_state.q2);
    dynamics.Jacobian_derivative[1][0] = -system_state.dq1 * sin(system_state.q1) - (system_state.dq1 + system_state.dq2) * sin(system_state.q1 + system_state.q2);
    dynamics.Jacobian_derivative[1][1] = -(system_state.dq1 + system_state.dq2) * sin(system_state.q1 + system_state.q2);

    dynamics.D_joint[0][0] = parameter.M[0] + parameter.M[1] + 2 * parameter.M[2] * cos(system_state.q2);
    dynamics.D_joint[0][1] = parameter.M[1] + parameter.M[2] * cos(system_state.q2);
    dynamics.D_joint[1][0] = parameter.M[1] + parameter.M[2] * cos(system_state.q2);
    dynamics.D_joint[1][1] = parameter.M[1];

    dynamics.C_joint[0][0] = -parameter.M[2] * system_state.dq2 * sin(system_state.q2);
    dynamics.C_joint[0][1] = -parameter.M[2] * (system_state.dq1 + system_state.dq2) * sin(system_state.q2);
    dynamics.C_joint[1][0] = parameter.M[2] * system_state.dq1 * sin(system_state.q2);
    dynamics.C_joint[1][1] = 0;

    dynamics.G_joint[0] = parameter.M[3] * parameter.g * cos(system_state.q1) + parameter.M[4] * parameter.g * cos(system_state.q1 + system_state.q2);
    dynamics.G_joint[1] = parameter.M[4] * parameter.g * cos(system_state.q1 + system_state.q2);

    double inv_Jacobi[OUT][OUT], inv_Jacobi_T[OUT][OUT], inv_Jacobi_T_D[OUT][OUT], inv_Jacobi_T_D_Jacobi_T[OUT][OUT], D_inv_Jacobi[OUT][OUT];
    double D_inv_Jacobi_dJ[OUT][OUT], C_D_inv_Jacobi_dJ[OUT][OUT], inv_Jacobi_T_C_D_inv_Jacobi_dJ[OUT][OUT], inv_D_task[OUT][OUT];
    double C_task_dx[OUT], Fx_C_task_dx_G_task[OUT], Jacobi_T[OUT][OUT];
    inv_matrix(inv_Jacobi, dynamics.Jacobian, 2);

    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            inv_Jacobi_T[k][j] = inv_Jacobi[j][k];
        }
    }

    matrix_Multi(inv_Jacobi_T_D, inv_Jacobi_T, dynamics.D_joint, OUT, OUT, OUT);
    matrix_Multi(dynamics.D_task, inv_Jacobi_T_D, inv_Jacobi, OUT, OUT, OUT);

    matrix_Multi(D_inv_Jacobi, dynamics.D_joint, inv_Jacobi, OUT, OUT, OUT);
    matrix_Multi(D_inv_Jacobi_dJ, D_inv_Jacobi, dynamics.Jacobian_derivative, OUT, OUT, OUT);

    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            C_D_inv_Jacobi_dJ[j][k] = dynamics.C_joint[j][k] - D_inv_Jacobi_dJ[j][k];
        }
    }

    matrix_Multi(inv_Jacobi_T_C_D_inv_Jacobi_dJ, inv_Jacobi_T, C_D_inv_Jacobi_dJ, OUT, OUT, OUT);
    matrix_Multi(dynamics.C_task, inv_Jacobi_T_C_D_inv_Jacobi_dJ, inv_Jacobi, OUT, OUT, OUT);

    matrix_Multi1(dynamics.G_task, inv_Jacobi_T, dynamics.G_joint, OUT, OUT, 1);

    inv_matrix(inv_D_task, dynamics.D_task, 2);

    for (int j = 0; j < OUT; j++){
        C_task_dx[j] = dynamics.C_joint[j][0] * system_state.dx1 + dynamics.C_joint[j][1] * system_state.dx2;
    }

    Fx_C_task_dx_G_task[0] = control_variable.force1 - C_task_dx[0] - dynamics.G_task[0];
    Fx_C_task_dx_G_task[1] = control_variable.force2 - C_task_dx[1] - dynamics.G_task[1];

    system_state.ddx1 = inv_D_task[0][0] * Fx_C_task_dx_G_task[0] + inv_D_task[0][1] * Fx_C_task_dx_G_task[1]; // dynamic equation of manipulator, Eq. 3.70
    system_state.ddx2 = inv_D_task[1][0] * Fx_C_task_dx_G_task[0] + inv_D_task[1][1] * Fx_C_task_dx_G_task[1];

    archive.ddx1_archive[i] = system_state.ddx1;
    archive.ddx2_archive[i] = system_state.ddx2;

    system_state.dx1 = system_state.dx1 + system_state.ddx1 * Ts;
    system_state.dx2 = system_state.dx2 + system_state.ddx2 * Ts;
    system_state.x1 = system_state.x1 + system_state.dx1 * Ts;
    system_state.x2 = system_state.x2 + system_state.dx2 * Ts;

    double sum1 = 0.0;
    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            sum1 += pow(dynamics.D_task[j][k], 2);
        }
    }
    dynamics.D_task_norm = sqrt(sum1);

    double sum2 = 0.0;
    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            sum2 += pow(dynamics.C_task[j][k], 2);
        }
    }
    dynamics.C_task_norm = sqrt(sum2);

    double sum3 = 0.0;
    for (int j = 0; j < OUT; j++){
        sum3 += pow(dynamics.G_task[j], 2);
    }
    dynamics.G_task_norm = sqrt(sum3);

    plant.plant_out1 = system_state.x1;
    plant.plant_out2 = system_state.dx1;
    plant.plant_out3 = system_state.x2;
    plant.plant_out4 = system_state.dx2;
    plant.plant_out5 = dynamics.D_task_norm;
    plant.plant_out6 = dynamics.C_task_norm;
    plant.plant_out7 = dynamics.G_task_norm;

    archive.D_task_norm_archive[i] = plant.plant_out5;
    archive.C_task_norm_archive[i] = plant.plant_out6;
    archive.G_task_norm_archive[i] = plant.plant_out7;

    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            Jacobi_T[k][j] = dynamics.Jacobian[j][k];
        }
    }

    control_variable.tol1 = Jacobi_T[0][0] * control_variable.force1 + Jacobi_T[0][1] * control_variable.force2;
    control_variable.tol2 = Jacobi_T[1][0] * control_variable.force1 + Jacobi_T[1][1] * control_variable.force2;

    archive.tol1_archive[i] = control_variable.tol1;
    archive.tol2_archive[i] = control_variable.tol2;
}

double CONTROL_realize(int i){
    controller.controller_u1 = xd1.y[i];
    controller.controller_u2 = dxd1.y[i];
    controller.controller_u3 = ddxd1.y[i];
    controller.controller_u4 = xd2.y[i];
    controller.controller_u5 = dxd2.y[i];
    controller.controller_u6 = ddxd2.y[i];
    controller.controller_u7 = system_state.x1;
    controller.controller_u8 = system_state.dx1;
    controller.controller_u9 = system_state.x2;
    controller.controller_u10 = system_state.dx2;
    // printf("system_state.x1 = %f\n", system_state.x1);
    archive.x1_archive[i] = controller.controller_u7;
    archive.dx1_archive[i] = controller.controller_u8;
    archive.x2_archive[i] = controller.controller_u9;
    archive.dx2_archive[i] = controller.controller_u10;

    for (int j = 0; j < H; j++){
        phi_D11[j] = exp((-pow(controller.controller_u7 - center_D[0][j], 2) - pow(controller.controller_u9 - center_D[1][j], 2)) / (width * width)); // output of the inertia matrix RBF function
        phi_D12[j] = exp((-pow(controller.controller_u7 - center_D[0][j], 2) - pow(controller.controller_u9 - center_D[1][j], 2)) / (width * width));
        phi_D21[j] = exp((-pow(controller.controller_u7 - center_D[0][j], 2) - pow(controller.controller_u9 - center_D[1][j], 2)) / (width * width));
        phi_D22[j] = exp((-pow(controller.controller_u7 - center_D[0][j], 2) - pow(controller.controller_u9 - center_D[1][j], 2)) / (width * width));
    }
    for (int j = 0; j < H; j++){
        phi_G1[j] = exp((-pow(controller.controller_u7 - center_G[0][j], 2) - pow(controller.controller_u9 - center_G[1][j], 2)) / (width * width)); // output of the gravity matrix RBF function
        phi_G2[j] = exp((-pow(controller.controller_u7 - center_G[0][j], 2) - pow(controller.controller_u9 - center_G[1][j], 2)) / (width * width));
    }
    for (int j = 0; j < H; j++){
        phi_C11[j] = exp((-pow(controller.controller_u7 - center_C[0][j], 2) - pow(controller.controller_u9 - center_C[1][j], 2) - pow(controller.controller_u8 - center_C[2][j], 2) - pow(controller.controller_u10 - center_C[3][j], 2)) / (width * width)); // output of the Coriolis matrix RBF function
        phi_C12[j] = exp((-pow(controller.controller_u7 - center_C[0][j], 2) - pow(controller.controller_u9 - center_C[1][j], 2) - pow(controller.controller_u8 - center_C[2][j], 2) - pow(controller.controller_u10 - center_C[3][j], 2)) / (width * width));
        phi_C21[j] = exp((-pow(controller.controller_u7 - center_C[0][j], 2) - pow(controller.controller_u9 - center_C[1][j], 2) - pow(controller.controller_u8 - center_C[2][j], 2) - pow(controller.controller_u10 - center_C[3][j], 2)) / (width * width));
        phi_C22[j] = exp((-pow(controller.controller_u7 - center_C[0][j], 2) - pow(controller.controller_u9 - center_C[1][j], 2) - pow(controller.controller_u8 - center_C[2][j], 2) - pow(controller.controller_u10 - center_C[3][j], 2)) / (width * width));
        // printf("phi_C11[%d] = %f\n", j, phi_C11[j]);
    }

    for (int j = 0; j < H; j++){
        weight_D11[j] = 0.0; // inertia matrix RBF network weight
        weight_D12[j] = 0.0;
        weight_D21[j] = 0.0;
        weight_D22[j] = 0.0;
        weight_G1[j] = 0.0; // gravity matrix RBF network weight
        weight_G2[j] = 0.0;
        weight_C11[j] = 0.0; // Coriolis matrix RBF network weight
        weight_C12[j] = 0.0;
        weight_C21[j] = 0.0;
        weight_C22[j] = 0.0;
    }

    controller.err1 = xd1.y[i] - controller.controller_u7;
    controller.err1_velocity = dxd1.y[i] - controller.controller_u8;
    controller.err2 = xd2.y[i] - controller.controller_u9;
    controller.err2_velocity = dxd2.y[i] - controller.controller_u10;
    archive.error1_archive[i] = controller.err1;
    archive.error1_velocity_archive[i] = controller.err1_velocity;
    archive.error2_archive[i] = controller.err2;
    archive.error2_velocity_archive[i] = controller.err2_velocity;
    // printf("controller.err1 = %f\n", controller.err1);
    controller.r1 = controller.err1_velocity + controller.Lambda[0][0] * controller.err1;
    controller.dxr1 = dxd1.y[i] + controller.Lambda[0][0] * controller.err1;
    controller.ddxr1 = ddxd1.y[i] + controller.Lambda[0][0] * controller.err1_velocity;
    controller.r2 = controller.err2_velocity + controller.Lambda[1][1] * controller.err2;
    controller.dxr2 = dxd2.y[i] + controller.Lambda[1][1] * controller.err2;
    controller.ddxr2 = ddxd2.y[i] + controller.Lambda[1][1] * controller.err2_velocity;
    // printf("controller.r1 = %f\n", controller.r1);

    // adaptive law
    for (int j = 0; j < H; j++){
        derivative_weight_D11[j] = controller.Tau_D11[j][j] * phi_D11[j] * controller.ddxr1 * controller.r1; // Eq. 3.60 and 3.77, adaptive law design
        derivative_weight_D12[j] = controller.Tau_D12[j][j] * phi_D12[j] * controller.ddxr2 * controller.r1;
        derivative_weight_D21[j] = controller.Tau_D21[j][j] * phi_D21[j] * controller.ddxr1 * controller.r2;
        derivative_weight_D22[j] = controller.Tau_D22[j][j] * phi_D22[j] * controller.ddxr2 * controller.r2;
    }
    for (int j = 0; j < H; j++){
        derivative_weight_G1[j] = controller.Tau_G1[j][j] * phi_G1[j] * controller.r1; // Eq. 3.62 and 3.77, adaptive law design
        derivative_weight_G2[j] = controller.Tau_G2[j][j] * phi_G2[j] * controller.r2;
    }
    for (int j = 0; j < H; j++){
        derivative_weight_C11[j] = controller.Tau_C11[j][j] * phi_C11[j] * controller.dxr1 * controller.r1; // Eq. 3.61 and 3.77, adaptive law design
        derivative_weight_C12[j] = controller.Tau_C12[j][j] * phi_C12[j] * controller.ddxr2 * controller.r1;
        derivative_weight_C21[j] = controller.Tau_C21[j][j] * phi_C21[j] * controller.dxr1 * controller.r2;
        derivative_weight_C22[j] = controller.Tau_C22[j][j] * phi_C22[j] * controller.ddxr2 * controller.r2;
        // printf("derivative_weight_C11[%d] = %f\n", j, derivative_weight_C11[j]);
    }

    for (int j = 0; j < H; j++){
        weight_D11[j] = weight_D11[j] + derivative_weight_D11[j] * Ts;
        weight_D12[j] = weight_D12[j] + derivative_weight_D12[j] * Ts;
        weight_D21[j] = weight_D21[j] + derivative_weight_D21[j] * Ts;
        weight_D11[j] = weight_D11[j] + derivative_weight_D11[j] * Ts;
    }
    for (int j = 0; j < H; j++){
        weight_G1[j] = weight_G1[j] + derivative_weight_G1[j] * Ts;
        weight_G2[j] = weight_G2[j] + derivative_weight_G2[j] * Ts;
    }
    for (int j = 0; j < H; j++){
        weight_C11[j] = weight_C11[j] + derivative_weight_C11[j] * Ts;
        weight_C12[j] = weight_C12[j] + derivative_weight_C12[j] * Ts;
        weight_C21[j] = weight_C21[j] + derivative_weight_C21[j] * Ts;
        weight_C22[j] = weight_C22[j] + derivative_weight_C22[j] * Ts;
    }

    for (int j = 0; j < H; j++){
        dynamics.DSNN_estimate[0][0] = weight_D11[j] * phi_D11[j]; // Eq. 3.51
        dynamics.DSNN_estimate[0][1] = weight_D12[j] * phi_D12[j];
        dynamics.DSNN_estimate[1][0] = weight_D21[j] * phi_D21[j];
        dynamics.DSNN_estimate[1][1] = weight_D22[j] * phi_D22[j];
    }
    for (int j = 0; j < H; j++){
        dynamics.GSNN_estimate[0] = weight_G1[j] * phi_G1[j]; // Eq. 3.53
        dynamics.GSNN_estimate[1] = weight_G2[j] * phi_G2[j];
    }
    for (int j = 0; j < H; j++){
        dynamics.CDNN_estimate[0][0] = weight_C11[j] * phi_C11[j]; // Eq. 3.52
        dynamics.CDNN_estimate[0][1] = weight_C12[j] * phi_C12[j];
        dynamics.CDNN_estimate[1][0] = weight_C21[j] * phi_C21[j];
        dynamics.CDNN_estimate[1][1] = weight_C22[j] * phi_C22[j];
    }

    double sum1 = 0.0;
    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            sum1 += pow(dynamics.DSNN_estimate[j][k], 2);
        }
    }
    dynamics.D_estimate_norm = sqrt(sum1);

    double sum2 = 0.0;
    for (int j = 0; j < OUT; j++){
        sum2 += pow(dynamics.GSNN_estimate[j], 2);
    }
    dynamics.G_estimate_norm = sqrt(sum2);

    double sum3 = 0.0;
    for (int j = 0; j < OUT; j++){
        for (int k = 0; k < OUT; k++){
            sum3 += pow(dynamics.CDNN_estimate[j][k], 2);
        }
    }
    dynamics.C_estimate_norm = sqrt(sum3);

    control_variable.force1 = dynamics.DSNN_estimate[0][0] * controller.ddxr1 + dynamics.DSNN_estimate[0][1] * controller.ddxr2 + dynamics.CDNN_estimate[0][0] * controller.dxr1 + dynamics.CDNN_estimate[0][1] * controller.dxr2 + dynamics.GSNN_estimate[0] + controller.K1 * controller.r1 + controller.Ks * (controller.r1 > 0 ? 1 : (controller.r1 < 0 ? -1 : 0)); // Eq. 3.75, controller design
    control_variable.force2 = dynamics.DSNN_estimate[1][0] * controller.ddxr1 + dynamics.DSNN_estimate[1][1] * controller.ddxr2 + dynamics.CDNN_estimate[1][0] * controller.dxr1 + dynamics.CDNN_estimate[1][1] * controller.dxr2 + dynamics.GSNN_estimate[1] + controller.K2 * controller.r2 + controller.Ks * (controller.r2 > 0 ? 1 : (controller.r2 < 0 ? -1 : 0)); // Eq. 3.75, controller design

    controller.controller_out1 = control_variable.force1;
    controller.controller_out2 = control_variable.force2;
    controller.controller_out3 = dynamics.D_estimate_norm;
    controller.controller_out4 = dynamics.G_estimate_norm;
    controller.controller_out5 = dynamics.C_estimate_norm;

    archive.Fx1_archive[i] = controller.controller_out1;
    archive.Fx2_archive[i] = controller.controller_out2;
    archive.D_task_estimate_norm_archive[i] = controller.controller_out3;
    archive.G_task_estimate_norm_archive[i] = controller.controller_out4;
    archive.C_task_estimate_norm_archive[i] = controller.controller_out5;
}

void saveArchiveToTxt(double *archive, int size, const char *filename){

    FILE *file = fopen(filename, "w+");

    if (file == NULL){
        perror("Failed to open file");
        exit(1);
    }
    else{
        for (int i = 0; i < size; i++){
            fprintf(file, "%lf\n", archive[i]);
        }
        fclose(file);
        printf("Saved to file %s\n", filename);
    }
}

void saveArchive(){

    saveArchiveToTxt(xd1.y, ARRAY_SIZE, "../report/xd1.txt");
    saveArchiveToTxt(archive.x1_archive, ARRAY_SIZE, "../report/x1.txt");
    saveArchiveToTxt(xd2.y, ARRAY_SIZE, "../report/xd2.txt");
    saveArchiveToTxt(archive.x2_archive, ARRAY_SIZE, "../report/x2.txt");
    saveArchiveToTxt(dxd1.y, ARRAY_SIZE, "../report/dxd1.txt");
    saveArchiveToTxt(archive.dx1_archive, ARRAY_SIZE, "../report/dx1.txt");
    saveArchiveToTxt(dxd2.y, ARRAY_SIZE, "../report/dxd2.txt");
    saveArchiveToTxt(archive.dx2_archive, ARRAY_SIZE, "../report/dx2.txt");
    saveArchiveToTxt(ddxd1.y, ARRAY_SIZE, "../report/ddxd1.txt");
    saveArchiveToTxt(archive.ddx1_archive, ARRAY_SIZE, "../report/ddx1.txt");
    saveArchiveToTxt(ddxd2.y, ARRAY_SIZE, "../report/ddxd2.txt");
    saveArchiveToTxt(archive.ddx2_archive, ARRAY_SIZE, "../report/ddx2.txt");
    saveArchiveToTxt(archive.error1_archive, ARRAY_SIZE, "../report/error1.txt");
    saveArchiveToTxt(archive.error1_velocity_archive, ARRAY_SIZE, "../report/error1_velocity.txt");
    saveArchiveToTxt(archive.error2_archive, ARRAY_SIZE, "../report/error2.txt");
    saveArchiveToTxt(archive.error2_velocity_archive, ARRAY_SIZE, "../report/error2_velocity.txt");
    saveArchiveToTxt(archive.Fx1_archive, ARRAY_SIZE, "../report/Fx1.txt");
    saveArchiveToTxt(archive.Fx2_archive, ARRAY_SIZE, "../report/Fx2.txt");
    saveArchiveToTxt(archive.tol1_archive, ARRAY_SIZE, "../report/tol1.txt");
    saveArchiveToTxt(archive.tol2_archive, ARRAY_SIZE, "../report/tol2.txt");
    saveArchiveToTxt(archive.D_task_estimate_norm_archive, ARRAY_SIZE, "../report/D_task_estimate_norm.txt");
    saveArchiveToTxt(archive.D_task_norm_archive, ARRAY_SIZE, "../report/D_norm.txt");
    saveArchiveToTxt(archive.G_task_estimate_norm_archive, ARRAY_SIZE, "../report/G_task_estimate_norm.txt");
    saveArchiveToTxt(archive.G_task_norm_archive, ARRAY_SIZE, "../report/G_task_norm.txt");
    saveArchiveToTxt(archive.C_task_estimate_norm_archive, ARRAY_SIZE, "../report/C_task_estimate_norm.txt");
    saveArchiveToTxt(archive.C_task_norm_archive, ARRAY_SIZE, "../report/C_task_norm.txt");
}

int main(){

    Input(&xd1, &dxd1, &ddxd1, &xd2, &dxd2, &ddxd2, Ts, t0, t1); // system input
    CONTROL_init();                                              // initialize controller parameter
    PLANT_init();                                                // initialize plant parameter

    for (int i = 0; i < ARRAY_SIZE; i++){
        double time = i * Ts + t0;
        printf("time at step %d: %f\n", i, time);
        CONTROL_realize(i);
        PLANT_realize(i);
    }

    saveArchive();

    return 0;
}
