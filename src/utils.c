#include <stdio.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "RKF78.h"

static float k_j[d_par] = {0.12, 0.0502, 0.0637, 0.1347, 0.0902, 0.0546,
                           0.0767, 0.1121,0.0971, 0.0403};
static float t_i[n_par] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 
                           27.0, 33.0};
static float eta_kj_cur[m_par][d_par] = {
    {0.0036, 0.0098, 0.0061, 0.0009, 0.0003, 0.0108, 0.0045, 0.0021, 0.0096, 0.0125},
    {0.0063, 0.0082, 0.0062, 0.0062, 0.0083, 0.013, 0.0039, 0.0019, 0.0015, 0.005},
    {0.0129, 0.0018, 0.0116, 0.0021, 0.009, 0.0129, 0.0054, 0.0049, 0.0093, 0.0066},
    {0.0053, 0.0086, 0.0067, 0.0029, 0.0089, 0.0054, 0.0042, 0.0095, 0.0112, 0.0092} 
};

static float eta_kj_pal[m_par][d_par] = {
    {0.00612, 0.01666, 0.01037, 0.00153, 0.00051, 0.01836, 0.00765, 0.00357, 0.01632, 0.02125},
    {0.01071, 0.01394, 0.01054, 0.01054, 0.01411, 0.0221, 0.00663, 0.00323, 0.00255, 0.0085},
    {0.02193, 0.00306, 0.01972, 0.00357, 0.0153, 0.02193, 0.00918, 0.00833, 0.01581, 0.01122},
    {0.00901, 0.01462, 0.01139, 0.00493, 0.01513, 0.00918, 0.00714, 0.01615, 0.01904, 0.01564}
};

void Gompertz(double t, double N, double *der, void *Params) {
    *der = ((N < 1.e-16) ? 0.0 :
            N*(lambdalogTheta_par - lambda_par*log(N) 
               - ((ODE_Parameters *) Params)->drift_i));
}

double Curative_Fitness(unsigned char *Cij) {
    register unsigned char i, j;
    ODE_Parameters GompertzParams;
    double N = NZero_par, t = t_i[0];
    double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
    unsigned char curativecounter = 0U, npar = n_par - 1;
    double integral = 0.0, lastt = t, lastN = N;

    if ( !TestIfConstraints2and3AreVerifiedCurative(Cij) ) return MAXDOUBLE;

    for ( i = 0; i < npar; i++ ) { double tfin = t_i[i+1]; 
        GompertzParams.drift_i = 0.0;
        for ( j = 0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

        while (t+h < tfin) {
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            if (N > NMax_par) return MAXDOUBLE;
            integral += (lastN + N) * (t - lastt);
            lastt = t; lastN = N;
        }

        do { h = tfin - t;
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            if (N > NMax_par) return MAXDOUBLE;
            integral += (lastN + N) * (t - lastt);
            lastt = t; lastN = N;
        } while (t < tfin);
        if (N < 1000) { curativecounter++; if (curativecounter > 2) return integral/2.0; }
        else curativecounter = 0U;
    }
    return MAXDOUBLE;
}

void writeCurativeSolution(unsigned char *Cij) {
    register unsigned char i, j;
    ODE_Parameters GompertzParams;
    double N = NZero_par, t = t_i[0];
    double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
    unsigned char curativecounter = 0U, npar = n_par - 1;
    double integral = 0.0, lastt = t, lastN = N;
    
    /* Opening the solutions file */
    FILE *fp;
    if ( (fp = fopen("curative_solution.sol", "w")) == NULL) {
        ExitError("when opening file to write", 1);
    }

    /* The solution is curative */
    fprintf(fp, "# t N\n %lf %lf\n", t, N);
    for ( i = 0; i < npar; i++ ) { double tfin = t_i[i+1]; 
        GompertzParams.drift_i = 0.0;
        for ( j = 0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

        while (t+h < tfin) {
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            fprintf(fp, " %lf %lf\n", t, N);
            integral += (lastN + N) * (t - lastt);
            lastt = t; lastN = N;
        }

        do { h = tfin - t;
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            fprintf(fp, " %lf %lf\n", t, N);
            integral += (lastN + N) * (t - lastt);
            lastt = t; lastN = N;
        } while (t < tfin);
        if (N < 1000) { 
            curativecounter++; 
            if (curativecounter > 2) {
                /* Setting the rest of Cij to zero nd printing the*/
                /* concentrations */
                for (int k = (i + 1) * d_par; k < npar*d_par; k++) Cij[k] = 0;
                fprintf(fp, "# integral : %lf\n", integral/2.0);
                break;
             }
        }
        else curativecounter = 0U;
    }
    /* printing the concentrations Cij i, dose, j drug */
    for (i = 0; i < npar; i++) {
        fprintf(fp, "# Dose %d : ", i+1);
        for ( j = 0; j < d_par; j++) {
            fprintf(fp, " %2d ", Cij[i*d_par +j]);
        }
        fprintf(fp, "\n");
    }       
}

double Paliative_Fitness(unsigned char *Cij) {
    register unsigned char i, j;
    ODE_Parameters GompertzParams;
    double N = NZero_par, t = t_i[0];
    double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
    unsigned char npar = n_par - 1;

    if ( !TestIfConstraints2and3AreVerifiedPaliative(Cij) ) return MAXDOUBLE;
            
    for ( i = 0; i < npar; i++ ) { double tfin = t_i[i+1]; 
        GompertzParams.drift_i = 0.0;
        for ( j = 0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

        while (t+h < tfin) {
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            if (N > NMax_par) return MAXDOUBLE;
        }

        do { h = tfin - t;
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            if (N > NMax_par) return MAXDOUBLE;
        } while (t < tfin);
    }
    /* After treatment is finished, no dose or treatment neither at time */ 
    /* tau_n nor later on */
    h = 1.e-3;
    GompertzParams.drift_i = 0.0;
    do {
        RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
    } while (N < NMax_par);
    return 1./t;
}

void writePaliativeSolution(unsigned char* Cij) {
    register unsigned char i, j;
    ODE_Parameters GompertzParams;
    double N = NZero_par, t = t_i[0];
    double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
    unsigned char npar = n_par - 1;

    /* opening the paliative solution file */
    FILE* fp;
    if ( (fp = fopen("paliative_solution.sol", "w")) == NULL) {
        ExitError("when opening file to write", 1);
    }

    /* The solution is paliative */ 
    fprintf(fp, "# t N\n %lf %lf\n", t, N);
    for ( i = 0; i < npar; i++ ) { double tfin = t_i[i+1]; 
        GompertzParams.drift_i = 0.0;
        for ( j = 0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * *(Cij++);

        while (t+h < tfin) {
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            fprintf(fp, " %lf %lf\n", t, N);
        }

        do { h = tfin - t;
            RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
            fprintf(fp, " %lf %lf\n", t, N);
        } while (t < tfin);
    }
    /* After treatment is finished, no dose or treatment neither at time */ 
    /* tau_n nor later on */
    h = 1.e-3;
    GompertzParams.drift_i = 0.0;
    do {
        RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
        fprintf(fp, " %lf %lf\n", t, N);
    } while (N < NMax_par);
    /* printing the doses and concentrations */
    for (i = 0; i < npar; i++) {
        fprintf(fp, "# Dose %d : ", i+1);
        for ( j = 0; j < d_par; j++) {
            fprintf(fp, " %2d ", Cij[i*d_par +j]);
        }
        fprintf(fp, "\n");
    }       
}

unsigned char TestIfConstraints2and3AreVerifiedCurative(unsigned char *Cij) {
    register unsigned char i, j, k;
    unsigned char npar = n_par - 1, *Cijofi;
    
    for ( j = 0; j < d_par; j++) { unsigned int ccumj = 0U;
        for ( i = 0; i < npar; i++) ccumj += *(Cij + i*d_par + j);
        if(ccumj > CCUMj_par) return 0U;
    } 
    
    for ( i = 0, Cijofi = Cij; i < npar; i++, Cijofi += d_par) {
        for ( k = 0; k < m_par; k++) { double Cseffk = 0.0;
            for( j = 0; j < d_par; j++) Cseffk += eta_kj_cur[k][j] * Cijofi[j];
            if(Cseffk > 1.0) return 0U;
        }
    }
    return 1U;
}

unsigned char TestIfConstraints2and3AreVerifiedPaliative(unsigned char *Cij) {
    register unsigned char i, j, k;
    unsigned char npar = n_par - 1, *Cijofi;
    
    for ( j = 0; j < d_par; j++) { unsigned int ccumj = 0U;
        for ( i = 0; i < npar; i++) ccumj += *(Cij + i*d_par + j);
        if(ccumj > CCUMj_par) return 0U;
    } 
    
    for ( i = 0, Cijofi = Cij; i < npar; i++, Cijofi += d_par) {
        for ( k = 0; k < m_par; k++) { double Cseffk = 0.0;
            for( j = 0; j < d_par; j++) Cseffk += eta_kj_pal[k][j] * Cijofi[j];
            if(Cseffk > 1.0) return 0U;
        }
    }
    return 1U;
}
