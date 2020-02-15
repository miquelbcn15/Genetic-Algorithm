#ifndef _UTILSH_
#define _UTILSH_

/* Tumour parameters */
#define lambda_par 0.336
#define lambdalogTheta_par 9.284023094951992774680543277273261947643
#define NZero_par 20000
#define CCUMj_par 127
#define NMax_par 9.5e+11
#define m_par 4
#define d_par 10
#define n_par 11
#define MAXDOUBLE DBL_MAX

typedef struct { double drift_i; } ODE_Parameters;
typedef double (*fitness)(unsigned char *Cij);
void Gompertz(double t, double N, double *der, void *Params);
double Curative_Fitness(unsigned char *Cij);
double Paliative_Fitness(unsigned char *Cij);
unsigned char TestIfConstraints2and3AreVerifiedCurative(unsigned char *Cij);
unsigned char TestIfConstraints2and3AreVerifiedPaliative(unsigned char *Cij);

#endif
