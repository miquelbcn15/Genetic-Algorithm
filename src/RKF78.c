/* Runge-Kutta-Fehlberg 78 with adaptive stepsize version 1.1
 * implemented by Lluis Alseda on Jan 15, 2020
 * This version corrects the error:
 *     bij += 2; // bij = &b_{i+1,1} and (bij += 2) = &b_{i+1,3}
 * to
 *     bij += 3; // bij = &b_{i+1,1} and (bij += 3) = &b_{i+1,4} */
#include "RKF78.h"
#include <math.h>

/* The functions RKF78 and RKF78Sys included here are an implementation
 * of the Runge-Kutta-Fehlberg method of orders 7 and 8 using a total of
 * s=13 stages (and evaluations of the vector field).
 * The difference between both estimates x^*_{n+1} and x_{n+1}
 * (with local errors of order 7 and 8) is computed and the L^1
 * norm is obtained. This norm is divided by n (the number of equations).
 * The number obtained in this way is required to be less than a given
 * tolerance tol times (1+0.01*nor), where nor is the L^1 norm of the
 * point computed to order 8. If this requirement is satisfied,
 * the order 8 estimate is taken as the next point.
 * Otherwise, a suitable value of the step h is obtained, and the
 * computation is started again.
 * In any case, when the next point is computed, a prediction of
 * the step h, to be used in the next call of the function, is found.
 *
 * Parameters:
 *   t: On input: a pointer to the current value of the independent
 *                variable; i.e., *t is the time corresponding
 *                to the actual initial condition.
 *      Output:   a pointer to the new time value *t corresponding to
 *                the new initial condition.
 *   x  position; i.e. the current value of the dependent variable
 *                (same input/output remarks as t).
 *                In RKF78Sys it is a vector of dimension n.
 *   h On input: a pointer to the time step *h to be used
 *               (it can be modified by the routine according to
 *               the given threshold).
 *     Output:   a pointer to the time step *h to be used in the next
 *               call of this function.
 *   hmin  the minimum allowed value for the absolute value of *h,
 *   hmax  the maximum allowed value for the absolute value of *h,
 *   tol   threshold to control the integration error.
 *   ParmsStruct: a pointer to a structure that contains the parameters
 *                of the ODE or vector field. It can be set to NULL if
 *                the vector field has no parameters.
 *                The format and name of the structure is free and
 *                decided by the user.
 *                The parameters of the ODE or the vector field (i.e.,
 *                the contents of the structure) are provided by the user.
 * For RKF78 it is necessary to provide
 *   ODE: a pointer (the name) of a function that computes the ODE.
 *            This function MUST be defined with the following header:

               void (*ODE)(double, double, double *, void *)

  or, for example (names are not compulsory),

               void ODEname(double t, double x, double *f, void *Params)

 *             where the first parameter (t -- input) is the time
 *                corresponding to the actual initial condition;
 *             the second parameter (x -- input) is the current initial
 *                condition (at time t), x(t);
 *             the third parameter (f -- output) is the evaluation of
 *                the ODE at time t on the point x=x(t)
 *             the fouth parameter (Params -- input) must be set to
 *                 * NULL if the ODE has no parameters, or
 *                 * the pointer to a structure that contains the
 *                   parameters of the ODE. Its format and
 *                   name are free and decided by conveniences of the
 *                   user.
 *
 * For RKF78Sys it is necessary to provide
 *   ODE_Sys: a pointer (the name) of a function that computes the vector
 *            field. This function MUST be defined with the following
 *            header:

               void (*ODE_Sys)(double, double *, unsigned char, double *, void *)

  or, for example (names are not compulsory),

               void vectorfieldname(double t, double x[], unsigned char n, double f[], void *Params)

 *             where the first parameter (t -- input) is the time
 *                corresponding to the actual initial condition;
 *             the second parameter (x -- input) is the current initial
 *                condition (at time t), x(t);
 *             the third parameter (n -- input) is the dimension of the
 *                vector field (and, hence, of x);
 *             the fouth parameter (f -- output -- a vector of dimension n)
 *                is the evaluation of the vector field derivatives at
 *                time t on the point x=x(t)
 *             the fifth parameter (Params -- input) must be set to
 *                 * NULL if the vector field has no parameters, or
 *                 * the pointer to a structure that contains the
 *                   parameters of the vector field. Its format and
 *                   name are free and decided by conveniences of the
 *                   user.
 * Returned value:
 *   an estimate of the error produced in the actual integration step.
 *
 * Unlike the function RKF78, the function RKF78Sys must be initialised
 * before its usage by the following call to the initialisation function:

              InitializeRKF78Sys(dimension);

 * Initialisation stores the dimension n of the vector field and creates
 * the following auxiliary storage (see also below):
 *   * a vector of dimension n to be used as workspace to store
 *       x_n + h \sum_{j=1}^{s-1} b_{s,j} k_j and the prediction of the RK7
 *   * a vector of dimension n to be used as workspace to store the prediction of the RK8
 *   * a matrix of dimension s x n to store the intermediate vectors k_i
*/

void RKF782tfin(double *t, double *x, double tfin, void *ParmsStruct,
              void (*ODE)(double, double, double *, void *))
{  double h;
	h = MIN(tfin - *t - 1.e-6, 1.e-3);
	while(1) { h = (*t + h < tfin) ? h : tfin - *t;
	   RKF78(t, x, &h, 1.e-8, 0.1, 1.e-14, ParmsStruct, ODE);
   }
}
void RKF78Sys2tfin(double *t, double *x, double tfin, void *ParmsStruct,
              void (*ODE_Sys)(double, double *, unsigned char, double *, void *))
{  double h;
	h = MIN(tfin - *t - 1.e-6, 1.e-3);
	while(1) { h = (*t + h < tfin) ? h : tfin - *t;
	   RKF78Sys(t, x, &h, 1.e-8, 0.1, 1.e-14, ParmsStruct, ODE_Sys);
   }
}

/* The Adaptive Runge-Kutta methods consist in applying the
 * same EXPLICIT Runge-Kutta method twice with different consecutive
 * orders.
 * This is achieved by using different coefficients for the different
 * evaluations of the vector field for every order.
 * These coefficients are represented in an extended Butcher tableau
 * of the form:

      a_1=0 |
      a_2   | b_{2,1}
      a_3   | b_{3,1} b_{3,2}
       .    |     .      .   .
       .    |     .      .     .
       .    |     .      .       .
      a_s   | b_{s,1} b_{s,2} ... b_{s,s-1}
     -------|--------------------------------------
            |   c_1     c_2   ...  c_{s-1}    c_s
            |  c^*_1   c^*_2  ... c^*_{s-1}  c^*_s

 * where the number s is called the number of stages.
 * We then compute

         x^*_{n+1} = x_n + h \sum_{i=1}^s c^*_i k_i
   and
           x_{n+1} = x_n + h \sum_{i=1}^s c_i k_i
 * where
 *   x^*_{n+1} is the solution given by the lower order method,
 *     x_{n+1} is the solution given by the higher order one, and
 * the k i’s are defined by:

        k_1 = f(t_n, x_n)
        k_2 = f(t_n + a_2 h, x_n + h b_{2,1} k_1)
        k_3 = f(t_n + a_3 h, x_n + h (b_{3,1} k_1 + b_{3,2} k_2))
            .
            .
            .
        k_s = f(t_n + a_s h, x_n + h \sum_{j=1}^{s-1} b_{s,j} k_j)

 * where f(t,x) denotes the vector field.
 *
 * What follows is an implementation of the extended Butcher tableau
 * for the RKF78
 * (with a total of s=13 stages (s=11 for the lower order method)).
 * The elements b_{i,j} with i=2,3,...,s and j=1,2,...,i-1
 * (note that  the number of elements b_{i,j} is

      1 + 2 + 3 + .... + (s-1 = 12) = (13 * 12)/2 = 13*6 = 78

 * ) are stored in the vector beta sorted lexicographically with respect
 * to the subindices i (first) and j (second); i.e. the vector beta is

     beta[0], beta[1], beta[2], ..., beta[66], beta[67], ..., beta[77]
     b_{2,1}, b_{3,1}, b_{3,2}, ..., b_{13,1}, b_{13,2}, ..., b_{13,12}

   or, in general,

      beta[i*(i-3)/2 + j] = b_{i,j} with i=2,3,...,s and j=1,2,...,i-1

 * The above nice analitic formula for the intex of beta in function of
 * the indices of b_{i,j}, is completely useless because if the
 * computation of the k_i's and the use of the k_j's in these
 * computations is done in increasing ordering of the i's (first) and
 * j's (second), then the appropriate successive indices of beta can be
 * computed simply by incrementing it by 1 starting from 0.
 *
 * NOTE on alpha: alpha[i] = a_{i-1} for i=0,1,2,...,12
 *    because the vectors in C start at 0.
 *    Moreover, since in the above formulae, alpha[0] is not used
 *    (it is included in the vector alpha to keep the aestetical fact
 *    that the element of the vector alpha corresponding to k_i is
 *    precisely alpha[i] (the subindices of k_i will also be renumbered
 *    from 0 to 12 to adapt the code to the C standard consisting on
 *    vectors and matrices starting at 0).
*/

static double alpha[13]={ 0.e0,
        2.e0/27.e0,      1.e0/9.e0,       1.e0/6.e0,      5.e0/12.e0,
             0.5e0,      5.e0/6.e0,       1.e0/6.e0,       2.e0/3.e0,
         1.e0/3.e0,           1.e0,            0.e0,            1.e0 };


static double beta[78]={
       2.e0/27.e0,
       1.e0/36.e0, 1.e0/12.e0,
       1.e0/24.e0,       0.e0,    1.e0/8.e0,
       5.e0/12.e0,       0.e0, -25.e0/16.e0,    25.e0/16.e0,
           0.5e-1,       0.e0,         0.e0,         0.25e0,           0.2e0,
    -25.e0/108.e0,       0.e0,         0.e0,  125.e0/108.e0,    -65.e0/27.e0,  125.e0/54.e0,
     31.e0/300.e0,       0.e0,         0.e0,           0.e0,    61.e0/225.e0,    -2.e0/9.e0,    13.e0/900.e0,
             2.e0,       0.e0,         0.e0,    -53.e0/6.e0,    704.e0/45.e0,  -107.e0/9.e0,     67.e0/90.e0,        3.e0,
    -91.e0/108.e0,       0.e0,         0.e0,   23.e0/108.e0,  -976.e0/135.e0,  311.e0/54.e0,    -19.e0/60.e0,  17.e0/6.e0,  -1.e0/12.e0,
  2383.e0/4100.e0,       0.e0,         0.e0, -341.e0/164.e0, 4496.e0/1025.e0, -301.e0/82.e0, 2133.e0/4100.e0, 45.e0/82.e0, 45.e0/164.e0, 18.e0/41.e0,
      3.e0/205.e0,       0.e0,         0.e0,           0.e0,            0.e0,   -6.e0/41.e0,    -3.e0/205.e0, -3.e0/41.e0,   3.e0/41.e0,  6.e0/41.e0, 0.e0,
 -1777.e0/4100.e0,       0.e0,         0.e0, -341.e0/164.e0, 4496.e0/1025.e0, -289.e0/82.e0, 2193.e0/4100.e0, 51.e0/82.e0, 33.e0/164.e0, 12.e0/41.e0, 0.e0, 1.e0 };

static double c7[11]={ 41.e0/840.e0,
	          0.e0,            0.e0,             0.e0,            0.e0,
	  34.e0/105.e0,      9.e0/35.e0,       9.e0/35.e0,     9.e0/280.e0,
	   9.e0/280.e0,    41.e0/840.e0 };

static double c8[13]={         0.e0,
	          0.e0,            0.e0,             0.e0,            0.e0,
	  34.e0/105.e0,      9.e0/35.e0,       9.e0/35.e0,     9.e0/280.e0,
       9.e0/280.e0,            0.e0,
      41.e0/840.e0,    41.e0/840.e0 };

/* Prototypes of utility functions */
double eighthroot(double);
#define MAX(a,b)	((a) > (b) ? (a) : (b))
#define ABS(x)		((x) < 0.0 ? -(x) : (x))

double veck[13];
double RKF78(double *t, double *x,
             double *h, double hmin, double hmax,
             double tol,
             void *ParmsStruct,
             void (*ODE)(double, double, double *, void *))
{	register unsigned char i, j;
	double err, tolr;
	double x7pred, x8pred;

/* Computing the values h*k_i with i=0,1,2,...,12
 * NOTE: The subindices have been numbered from 0 to 12
 *       (instead of from 1 to s=13) to adapt the code to the
 *       C standard consisting on vectors and matrices starting at 0
 * Taking this convention in to account and the notes to the
 * implementation of the extended Butcher tableau ve have:

     k_i = f(t_n + alpha[i]*h, x_n + \sum_{j=0}^{i-1} *(bij++) (h*k_j)) for i \geq 1

 * after initializing bij=beta.
*/
	while (1) { register double *bij;
/* Computation of h * k_0. Special case since x is not displaced */
		 ODE(*t, *x, veck, ParmsStruct); veck[0] *= *h; // vec[0] = k_0; veck = &vec[0] = &k_0

/* Computation of h * k_1 */
		 x7pred = *x + beta[0] * veck[0];
		 ODE(*t + alpha[1]* *h, x7pred, veck+1, ParmsStruct); veck[1] *= *h; // vec[i] = k_i; veck+i = &vec[i] = &k_i

/* Computation of h * k_2 */
         x7pred = *x + beta[1] * veck[0] + beta[2] * veck[1];
		 ODE(*t + alpha[2]* *h, x7pred, veck+2, ParmsStruct); veck[2] *= *h;

/* Computation of h * k_3. Special case: b_{4,2} = beta[4] = 0. To avoid multiplying by 0 */
         x7pred = *x + beta[3] * veck[0] + beta[5] * veck[2];
		 ODE(*t + alpha[3]* *h, x7pred, veck+3, ParmsStruct); veck[3] *= *h;

/* Computation of h * k_4. Special case: b_{5,2} = beta[7] = 0. To avoid multiplying by 0 */
         x7pred = *x + beta[6] * veck[0] + beta[8] * veck[2] + beta[9] * veck[3];
		 ODE(*t + alpha[4]* *h, x7pred, veck+4, ParmsStruct); veck[4] *= *h;

/* Computation of h * k_i; i >= 5. Normal cases: b_{i+1,2} = b_{i+1,3} = 0 */
		 bij=&beta[10];
		 for(i=5; i < 13; i++){
			 x7pred = *x + *bij * veck[0]; bij += 3; // bij = &b_{i+1,1} and (bij += 3) = &b_{i+1,4}
			 for(j=3; j < i; j++) x7pred += *(bij++) * veck[j];
			 ODE(*t + alpha[i]* *h, x7pred, veck+i, ParmsStruct); veck[i] *= *h;
		 }

/*    Computing the rk7 and rk8 predictions. We use that
 *                c7[1] = c7[2] = c7[3] = c7[4] = 0, and
 *        c8[0] = c8[1] = c8[2] = c8[3] = c8[4] = c8[10] = 0
*/
         x7pred = *x + c7[0] * veck[0]; x8pred = *x;
         for(j=5; j < 10; j++){ x7pred += c7[j] * veck[j]; x8pred += c8[j] * veck[j]; }
         x7pred += c7[10] * veck[10];
         x8pred += c8[11] * veck[11] + c8[12] * veck[12];

		 err = fabs(x8pred-x7pred); // Prediction error
         tolr = tol * (1.0 + fabs(x8pred)/100.0); /* It is an absolute
  tolerance for small values of the components and
  relative but "retarded two digts" for large values */

         if( ABS(*h) <= hmin || err < tolr) break; // Stop if we are already performing at minimum stepsize or the error is OK

         *h *= 0.9 * eighthroot(tolr/err); // Otherwise we recompute the step to try again
                             // NOTE that since err >= tolr ==> tolr/err <= 1.0 ==> 0.9 * pow(tolr/err, 0.125) <= 0.9 ==> fabs(new_h) < fabs(old_h) <= hmax
         *h = (*h < 0.0) ? (*h > - hmin ? -hmin : *h) : (*h < hmin ? hmin : *h); //Returning h to valid region, if necessary. Recall that fabs(*h) < hmax
    } // End of main loop; Going to recompute the RK78 estimates

     *t += *h; *x = x8pred; // Storing final time and x_{n+1}

 /* Step correction */
     err = MAX(err, tolr/256); // If the error was very small we do not want that the new step becomes too large (err >= tolr/256 ==> tolr/err <= 256)
     *h *= 0.9 * eighthroot(tolr/err); // Fehlberg correction (Stoer (7.2.5.16))
     *h = (*h < 0.0) ? (*h > - hmin ? -hmin : (*h < -hmax ? -hmax : *h)) : (*h < hmin ? hmin : (*h > hmax ? hmax : *h)); // Returning h to valid region, if necessary

    return (err);
} /* and this is the end */


/* BLAS-like convenience functions */
void VectorCopy_double(double *, short, double *);
void VectorAddCntntTimesVector_double(double *, short, double, double *);
void VectorSetToVectorPlusCntntTimesVector_double(double *, short, double *, double, double *);
void VectorSetToCntntTimesVector_double(double *, short, double, double *);
void VectorMultiplyByCntnt_double(double *, short, double);

double *matki, *x7pred, *x8pred; // matki[13][n]; x7pred[n]; x8pred[n]
unsigned char n;
unsigned short n2, n3;

double RKF78Sys(double *t, double x[],
             double *h, double hmin, double hmax,
             double tol,
             void *ParmsStruct,
             void (*ODE_Sys)(double, double *, unsigned char, double *, void *))
{	register unsigned char i, j;
	double err, nor, tolr;

/* Computing the vectors h*k_i with i=0,1,2,...,12
 * NOTE: The subindices have been numbered from 0 to 12
 *       (instead of from 1 to s=13) to adapt the code to the
 *       C standard consisting on vectors and matrices starting at 0
 * Taking this convention in to account and the notes to the
 * implementation of the extended Butcher tableau ve have:

     k_i = f(t_n + alpha[i]*h, x_n + \sum_{j=0}^{i-1} *(bij++) (h*k_j)) for i \geq 1

 * after initializing bij=beta.
 *
 * NOTE: beta[m++] = *(beta++) by using pointer arithmetic.
 * NOTE: h*k_i stored in the vector &matki[i][0]
 * Using the vector x7pred as intermediate storage
*/
	while (1) { register double *kj, *bij;
/* Computation of h * k_0. Special case since x is not displaced */
		 ODE_Sys(*t, x, n, matki, ParmsStruct); // k_0 = f(t, x) k_0 -> &matki[i][0] = matki
		 VectorMultiplyByCntnt_double(matki, n, *h); // k_0 -> h * k_0

/* Computation of h * k_1 */
		 VectorSetToVectorPlusCntntTimesVector_double(x7pred, n, x, beta[0], matki); // x7pred = x + beta[0] (h * k_0)
		 kj=matki+n; // k_1 -> &matki[1][0] = matki+n
		 ODE_Sys(*t + alpha[1]* *h, x7pred, n, kj, ParmsStruct); // k_1 = f(t + alpha[1]*h, x + beta[0]*(h * k_0));
		 VectorMultiplyByCntnt_double(kj, n, *h); // k_1 -> h * k_1

/* Computation of h * k_2 */
		 VectorSetToVectorPlusCntntTimesVector_double(x7pred, n, x, beta[1], matki); // x7pred = x + beta[1] (h * k_0)
		 VectorAddCntntTimesVector_double(x7pred, n, beta[2], kj); // x7pred += beta[2]*(h * k_1) since kj=matki+n
		 kj=matki+n2; // k_2 -> &matki[2][0] = matki+2*n
		 ODE_Sys(*t + alpha[2]* *h, x7pred, n, kj, ParmsStruct); // k_2 = f(t + alpha[2]*h, x + beta[1] (h * k_0) + beta[2] (h * k_1));
		 VectorMultiplyByCntnt_double(kj, n, *h); // k_2 -> h * k_2

/* Computation of h * k_3. Special case: b_{4,2} = beta[4] = 0. To avoid multiplying by 0 */
		 VectorSetToVectorPlusCntntTimesVector_double(x7pred, n, x, beta[3], matki); // x7pred = x + beta[3] (h * k_0)
		 VectorAddCntntTimesVector_double(x7pred, n, beta[5], kj); // x7pred += beta[5]*(h * k_2) since kj=matki+2*n
		 kj=matki+n3; // k_3 -> &matki[3][0] = matki+3*n
		 ODE_Sys(*t + alpha[3]* *h, x7pred, n, kj, ParmsStruct); // k_3 = f(t + alpha[3]*h, x + beta[3] (h * k_0) + beta[5] (h * k_2));
		 VectorMultiplyByCntnt_double(kj, n, *h); // k_3 -> h * k_3

/* Computation of h * k_4. Special case: b_{5,2} = beta[7] = 0. To avoid multiplying by 0 */
		 VectorSetToVectorPlusCntntTimesVector_double(x7pred, n, x, beta[6], matki); // x7pred = x + beta[6] (h * k_0)
		 VectorAddCntntTimesVector_double(x7pred, n, beta[8], matki+n2); // x7pred += beta[8]*(h * k_2)
		 VectorAddCntntTimesVector_double(x7pred, n, beta[9], kj); // x7pred += beta[9]*(h * k_3) since kj=matki+3*n
		 kj += n; // k_4 -> &matki[4][0] = matki+4*n
		 ODE_Sys(*t + alpha[4]* *h, x7pred, n, kj, ParmsStruct); // k_3 = f(t + alpha[4]*h, x + beta[6] (h * k_0) + beta[8] (h * k_2) + beta[9]*(h * k_3));
		 VectorMultiplyByCntnt_double(kj, n, *h); // k_4 -> h * k_4

/* Computation of h * k_i; i >= 5. Normal cases: b_{i+1,2} = b_{i+1,3} = 0 */
		 bij=&beta[10];
		 for(i=5; i < 13; i++){
		 	 VectorSetToVectorPlusCntntTimesVector_double(x7pred, n, x, *bij, matki); // x7pred = x + *bij * (h*k_0); bij = &b_{i+1,1}
		 	 bij += 3; // bij = &b_{i+1,4}
			 for(j=3, kj=matki+n3; j < i; j++, kj+=n) VectorAddCntntTimesVector_double(x7pred, n, *(bij++), kj); // x7pred += *(bij++) * (h*k_j)
			 ODE_Sys(*t + alpha[i]* *h, x7pred, n, kj, ParmsStruct); // k_i -> &matki[i][0] = matki+i*n = kj
			 VectorMultiplyByCntnt_double(kj, n, *h); // k_i -> h * k_i
		 }

/*    Computing the rk7 and rk8 predictions. We use that
 *                c7[1] = c7[2] = c7[3] = c7[4] = 0, and
 *        c8[0] = c8[1] = c8[2] = c8[3] = c8[4] = c8[10] = 0
 *    For speed we priorize to run matki sequentially (incrementally by 1)
*/
         err = nor = 0.0;
		 for(i=0, kj = matki; i < n; i++) { x8pred[i] = x[i]; x7pred[i] = x[i] + *(kj++) * c7[0]; } // k_0 -> &matki[i][0] = matki; x*pred = x + c*[0]*(h*K_0)
	     for(j=5, kj = matki + 5*n; j < 10; j++){ for(i=0; i < n; i++, kj++) { x7pred[i] += *kj * c7[j]; x8pred[i] += *kj * c8[j]; } } // x*pred = += c*[j]*(h*K_j)
	     VectorAddCntntTimesVector_double(x7pred, n, c7[10], kj);
	     VectorAddCntntTimesVector_double(x8pred, n, c8[11], kj+n);
	     for(i=0, kj += n2; i < n; i++, kj++) { x8pred[i] += *kj * c8[12];
			 nor += fabs(x8pred[i]); // L^1 norm of the prediction. Used to relativize tolerance
			 err += fabs(x8pred[i]-x7pred[i]);
		 }

		 err /= n; // average of the component errors
         tolr = tol * (1.0 + nor/100.0); /* It is an absolute
  tolerance for small values of the components and
  relative but "retarded two digts" for large values */

         if( ABS(*h) <= hmin || err < tolr) break; // Stop if we are already performing at minimum stepsize or the error is OK

         *h *= 0.9 * eighthroot(tolr/err); // Otherwise we recompute the step to try again
                    // NOTE that since err >= tolr ==> tolr/err <= 1.0 ==> 0.9 * pow(tolr/err, 0.125) <= 0.9 ==> fabs(new_h) < fabs(old_h) <= hmax
         *h = (*h < 0.0) ? (*h > - hmin ? -hmin : *h) : (*h < hmin ? hmin : *h); //Returning h to valid region, if necessary. Recall that fabs(*h) < hmax
    } // End of main loop; Going to recompute the RK78 estimates

     *t += *h; // Storing final time
     VectorCopy_double(x, n, x8pred); // Storing x_{n+1}

 /* Step correction */
     err = MAX(err, tolr/256); // If the error was very small we do not want that the new step becomes too large (err >= tolr/256 ==> tolr/err <= 256)
     *h *= 0.9 * eighthroot(tolr/err); // Fehlberg correction (Stoer (7.2.5.16))
     *h = (*h < 0.0) ? (*h > - hmin ? -hmin : (*h < -hmax ? -hmax : *h)) : (*h < hmin ? hmin : (*h > hmax ? hmax : *h)); // Returning h to valid region, if necessary

    return (err);
} /* and this is the end */

/* Initialisation of RKF78
 * This function stores the dimension n of the vector field and creates
 * the following auxiliary storage:
 *   * a vector od dimennsion n to store x_n + h \sum_{j=1}^{s-1} b_{s,j} k_j
 *   * a matrix of dimension s x n to store the intermediate vectors k_i      */
void ExitError(const char *miss, int errcode) { fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss); exit(errcode); }
void InitializeRKF78Sys(unsigned char dimension){ n = dimension; n2 = n + n; n3 = n2 + n;
	if((matki = (double *) malloc(13*dimension*sizeof(double))) == NULL) ExitError("when allocating memory for the K_i's matrix", 3);
	if((x7pred = (double *) malloc(dimension*sizeof(double))) == NULL) ExitError("when allocating internal storage space", 5);
	if((x8pred = (double *) malloc(dimension*sizeof(double))) == NULL) ExitError("when allocating internal storage space", 5);
}

/* eighthrootofpowersoftwo[i] = 2^(i/8) */
static double eighthrootofpowersoftwo[] = { 1.0,
    1.090507732665257659207010655760707978993,
    1.189207115002721066717499970560475915293,
    1.296839554651009665933754117792451159836,
    1.41421356237309504880168872420969807857,
    1.542210825407940823612291862090734841307,
    1.68179283050742908606225095246642979008,
    1.834008086409342463487083189588288856078 };

static double twoto64 = 1.8446744073709551616e19;

/* scalerootfactor[i] = ((256+i)/256)^(1/8) */
static double scalerootfactor[] = { 1.0,
	1.000487448816538668711013867569244074793, 1.00097324084707034743604501683773495937,
    1.001457388111856589232056240986421378863, 1.001939902497953687584112101455088116314,
    1.002420795761194625568019539246451962849, 1.002900079528134067496301843512425335846,
    1.003377765297957220396743200994905202723, 1.003853864444353371080950508389247131301,
    1.00432838821735488361160406240505651945,  1.00480134774514242164980940407031401131,
    1.005272754035817140437516098405269381273, 1.005742617979140574024397146865066186677,
    1.006210950348242924764680111670836216809, 1.006677761801300444068686811666383312895,
    1.00714306288318257587844713631089644477,  1.007606864027069517329527132872536563458,
    1.00806917555604083454559768708206382991,  1.008530007684635755472311332060034890454,
    1.008989370520385746077368642869627757241, 1.009447274065319961109411459268201301031,
    1.00990372821744414590527678183107918786,  1.010358742772193551449390593879966363477,
    1.010812327423860411007372211728970888127, 1.01126449176699651316542431091966587886,
    1.011715245297791392995420652028605006481, 1.012164597415426650320826293956135592606,
    1.012612557423406891669164987725584034197, 1.013059134530867780451558410398989725184,
    1.013504337853861668198161158350767644044, 1.013948176416621268289734795950937043811,
    1.014390659152801822550132058103191045406, 1.014831794906702200292430800401869400744,
    1.015271592434465358933529736681669326622, 1.015710060405258585099176212193169210408,
    1.016147207402433925224928669356403213492, 1.016583041924669205010046693972076897133,
    1.017017572387090027692617411705594514757, 1.017450807122373131977510040210230128422,
    1.017882754381831481556405659608448689025, 1.018313422336481449503835675042560344248,
    1.018742819078092452407783491856357883342, 1.019170952620219380891098697837976595717,
    1.019597830899218165194107662292016488159, 1.020023461775244806712963811100259094355,
    1.020447853033238198816260671834525625722, 1.020871012383887052888230106458874531313,
    1.021292947464581238364661907275508427049, 1.02171366584034783853189299106603545184,
    1.022133175004772217044390825229099123333, 1.0225514823809043834773381820781865576,
    1.022968595322150939762126018001162280289, 1.023384521113152883049852844198615325319,
    1.023799266970649534406044719874111515872, 1.02421284004432885675423443296185206015,
    1.024625247417664419652302834810707551527, 1.025036496108739262799262655044472041447,
    1.025446593071056904627265194709867105583, 1.025855545194339736929974850560083833398,
    1.026263359305315041210154682846649959611, 1.0266700421684888572925304823928465027,
    1.027075600486907929739062127077274655277, 1.027480040902909952719075307299006694774,
    1.027883369998862329222830646584928005691, 1.028285594297889655860674613788187804615,
    1.028686720264590139957674614768590444409, 1.029086754305741151232436393044085955049,
    1.029485702770994106035579050733739054849, 1.029883571953558877915138627600505390051,
    1.030280368090877924170112438029143320849, 1.030676097365290314046657626928507742045,
    1.03107076590468584032141729800303252527,  1.031464379783149392200446193561616590844,
    1.031856945021595763737704138180219088073, 1.032248467588395068341614345093453964899,
    1.032638953399988926389353935065461347347, 1.033028408321497589504035591415383132773,
    1.033416838167318161667501047662783399774, 1.033804248701714074038894621068837699853,
    1.034190645639395967124398312748197170973, 1.034576034646094130794431549815678587528,
    1.034960421339122649569251296878853017277, 1.035343811287935397590293309253538405005,
    1.035726210014674024760890603141666367564, 1.03610762299470807267436339929008766579,
    1.036488055657167356148121483852802443257, 1.036867513385466743447632069384903924546,
    1.037246001517823465612210464978784024667, 1.037623525347767082683961927765640253165,
    1.038000090124642232090262314809213402065, 1.038375701054104281937379060535378084078,
    1.038750363298608009536712758588930651272, 1.039124081977889423104235761059672701996,
    1.039496862169440842246611277787574845928, 1.039868708908979350572827765398388270225,
    1.040239627190908731546650753900660631355, 1.040609621968774996521486803006004979729,
    1.040978698155715611774117314976582727672, 1.04134686062490252927597381625366069501,
    1.041714114209979123909005397589101824332, 1.042080463705491137846580554125313195161,
    1.042445913867311730877148911065141219783, 1.042810469413060733548473439779633241035,
    1.043174135022518198152070970811839536754, 1.043536915338032340750037418069284157521,
    1.043898814964921965668681704305007171501, 1.044259838471873462144373871808552295848,
    1.044619990391332461105779790600535455213, 1.044979275219890238412284527433507579687,
    1.045337697418664949240001103140516191981, 1.045695261413677776713447592388480658415,
    1.046051971596224076321903396248510930147, 1.046407832323239596133797973594583241378,
    1.046762847917661851329437449215418194584, 1.047117022668786730111152936780906617376,
    1.047470360832620406619796619004728362405, 1.047822866632226635086675383067770406385,
    1.048174544258069498079774564271338825434, 1.048525397868351680361782659598724340283,
    1.048875431589348338564296842028299360099, 1.049224649515736635597001857350879095256,
    1.049573055710921007451922033394731643727, 1.049920654207354228830415282596193778887,
    1.05026744900685434281379324031924037931,  1.050613444080917518616713189712291214363,
    1.050958643371026900305210886426515263961, 1.051303050788957508227859647247001279216,
    1.051646670217077253798495636305552876495, 1.051989505508644127181701981581981765079,
    1.052331560488099616367268869269211011133, 1.05267283895135841507663026505839139667,
    1.053013344666094475922320675172490244832, 1.053353081372023464240310398771444681881,
    1.053692052781181667034190430352950812273, 1.054030262578201410509125973891835152567,
    1.054367714420583038731829557126652940806, 1.054704411938963505030081486520108338211,
    1.055040358737381626841118422166828381108, 1.055375558393540053832102492054604951215,
    1.055710014459063998247466378451917481502, 1.05604373045975677558680713715114867035,
    1.056376709895852202882785984736462859879, 1.056708956242263901030805361763503343149,
    1.057040472948831546820710059864546149058, 1.057371263440564119535037000028605726981,
    1.057701331117880186208068134716201676594, 1.058030679356845268884781300587171565988,
    1.058359311509406336478411436634309526832, 1.058687230903623463099404326303038083792,
    1.059014440843898694016749780976786857571, 1.059340944611202159714711539715618422751,
    1.059666745463295477823525217821825939601, 1.059991846634952482031418813254362657465,
    1.060316251338177316427035115995630967829, 1.060639962762419933075721335173731271255,
    1.0609629840747890299999245854561156189,   1.061285318420262466112825352363351685044,
    1.061606968921895189045093877859146765076, 1.061927938681024711207011994471251758742,
    1.062248230777474168841916775900728639541, 1.062567848269752998251749859402719845355,
    1.062886794195255262811200577929834071933, 1.06320507157045566383328086894475548838,
    1.063522683391103267805939510200491205029, 1.063839632632412981986292094867975608622,
    1.064155922249254809815996001816477703969, 1.064471555176340917108026189090938054464,
    1.064786534328410539451402602092382375264, 1.065100862600412760786082793775605204153,
    1.065414542867687192615068119469515678298, 1.065727577986142582844587257578416771544,
    1.066039970792433382775829914105359046212, 1.066351724104134300312923806379410743535,
    1.066662840719912867001501010748904118537, 1.066973323419700046070111223383185705975,
    1.067283174964858908212739144011322159419, 1.06759239809835140142460466818518610037,
    1.067900995544903240785105282701815302729, 1.068208970011166943671041123490020890558,
    1.068516324185883035479989309288282616401, 1.06882306074003945054771366274492247518,
    1.069129182327029152554660461547259555148, 1.06943469158280599833475546232504212613,
    1.069739591126038868624740410920543914119, 1.070043883558264088924030079970850251598,
    1.070347571464036163273398150032568390883, 1.07065065741107684340557959424006433618,
    1.070953143950422555371979213933616613596, 1.071255033616570205406974061234005161442,
    1.071556328927621386454667943262211397481, 1.071857032385425006452278043725827010675,
    1.072157146475718359139488609248440860448, 1.072456673668266657843978929602041723693,
    1.072755616417001052379809340255631499499, 1.073053977160155148887319025908749100168,
    1.073351758320400052140544764637122505613, 1.073648962304977949550804551361553025336,
    1.073945591505834255802900714792461475823, 1.074241648299748336773282387421557298854,
    1.074537135048462831097367898564016264372, 1.074832054098811587475966879783501684033,
    1.075126407782846235538264741793037296744, 1.075420198417961407811045891433137525214,
    1.075713428307018630080645795104538491584, 1.076006099738468897175446901478077778135,
    1.076298214986473950942482557600947028793, 1.076589776311026276941801295714356289141,
    1.07688078595806783613658795734410278946,  1.077171246159607547615556555137676495448,
    1.077461159133837538146742783027923936901, 1.077750527085248174128453595556079767447,
    1.078039352204741891273700864502801491383, 1.078327636669745837138880976779183690026,
    1.07861538264432334138568913027113640559,  1.07890259227928422844720432295300691851,
    1.079189267712293987054678428913984063746, 1.079475411067981810870741596964400438197,
    1.079761024458047524267429213524146593877, 1.080046109981367407083576965189519763196,
    1.080330669724098931995655612540148359788, 1.08061470575978442793896277999735229186,
    1.080898220149453682822193520711604690617, 1.081181214941725498587714053769494576667,
    1.081463692172908211482304567016087791114, 1.081745653867099190218659227988817455749,
    1.082027102036283324526477637651757618675, 1.082308038680430516413496153920274056121,
    1.082588465787592186281235202496154182795, 1.082868385333996805867526394750367772624,
    1.083147799284144469817978590051986623532, 1.083426709590900517521393642678771417499,
    1.083705118195588216679700172933451446885, 1.083983027028080519921188028801665501538,
    1.084260438006890905606648885007270360931, 1.08453735303926331382141236121199948342,
    1.084813774021261188392165783315186363649, 1.085089702837855635615813845588223842881,
    1.08536514136301271023842745172013770625,  1.085640091459779839075505305415263623274,
    1.085914554980371392520284642709263043301, 1.086188533766253414044646957350408152254,
    1.086462029648227517657229609417974586931, 1.086735044446513963145634582883832390955,
    1.087007579970833918794081924424847877388, 1.087279638020490921134448887752499496187,
    1.087551220384451541157328622315469408126, 1.087822328841425266280497226160113852059,
    1.088092965159943607244958703991477257582, 1.08836313109843843898350812155716033488,
    1.088632828405319584383479013971738254165, 1.088902058819051649743987560163078959209,
    1.089170824068230120608519519275813672038, 1.089439125871656726536093433657333847919,
    1.08970696593841408325844277403758721041,  1.089974345967939620556658800084936154018,
    1.090241267650098804078493812030024120755 };

/* Computes the eighth root of a number x \in (0, 256]
 * Algorithm:
 *  Stage 1: For any double normal number x \in (0, 256] there exist unique
 *     p \in \{256,257,258,...,511\}, z \in [0, 1/p) and expnt
 *     such that
          x = (p/256) * (1 + z) * 2^expnt =
            = (p/256) * (1 + z) * 2^(8*int(expnt/8)) * 2^(expnt%8)
 *     Moreover, -1022 <= expnt <= floor(log_2(256)) = 8
 *     Proof: Let expnt := floor(log_2(x)).
 *        Then log_2(x) = expnt + alpha with alpha \in [0,1). Thus,
 *        x = 2^(log_2(x)) = 2^expnt * 2^(alpha) with 2^(alpha) \in [1,2).
 *        Now, set
 *           p := floor(256*2^(alpha)) \in \{256,257,258,...,511\},  and
 *           z := 256*2^(alpha)/p - 1 (that is, 2^(alpha) = (p/256) *(1 + z))
 *        Clearly,   p <= 256*2^(alpha) < p+1 ==>
 *                   1 <= 256*2^(alpha)/p < 1 + 1/p ==>
 *                   0 <= z < 1/p
 *        p/256 <= 2^(alpha) < (p+1)/256.
 *        Moreover, expnt <= floor(log_2(256)) = 8.
 *        Also, since a double normal number is of the form
               x = 2^e * 1.fraction with e >= -1022
 *        we have -1022 <= e = floor(log_2(2^e * 1.fraction)) = expnt
 *                                                                Q.E.D.
 * Consequently, since 256 = 2^8,
     x^(1/8) = (p/256)^(1/8) * (1 + z)^(1/8) * 2^(int(expnt/8)) * 2^((expnt%8)/8)
 *
 *  Stage 2: Approximate (1+z)^(1/8) by
        1 + z/8  - (7*z^2)/128  + (35*z^3)/1024 - (805*z^4)/32768 + (4991*z^5)/262144 - (64883*z^6)/4194304 =
        1 + z*(1 + z*(z*(35 + z*(z*(4991 - 4055.1875*z)/8 - 805)/32)/8 - 7)/16)/8
 *  because 4194304 = 262144 * 16 = 32768 * 8 * 16 = 1024 * 32 * 8 * 16 = 128 * 8 * 32 * 8 * 16	= 8 * 16 * 8 * 32 * 8 * 16;
 *  and 64883/16 = 4055.1875
 *
 * Test 512e+6 evaluations of pow(x, 0.125): 0m24.847s
 *                            eighthroot(x): 0m0.669s
*/
double eighthroot(double x){
	union { double val;
            struct { unsigned long  int  mantisa : 52;
                     unsigned short int  exp  : 11; // From 0 to 2047 biased by 1023
                     unsigned       char sign : 1; } bits;
	} doublesplittable; doublesplittable.val = x;

	if(x < 0.0 || doublesplittable.bits.exp > 1038U) return (0.0/0.0); // In fact the function computes x^(1/8) for 2^e * 1.fraction = x with e <= 15 (that is, for x < 2^15 * 2 = 2^16
	if(!doublesplittable.bits.exp) return 0.0; // x is not a normal number

/* Extracting the exponent of two of x. Since x is normal is of the form
 * x = 2^(expnt) * 1.fraction; and this is the same as computing expnt = floor(log_2(x)) */
	int expnt = ((int) doublesplittable.bits.exp) - 1023;
	doublesplittable.bits.exp = 1023U; // Setting the exponent of two of x to zero so doublesplittable.val is now 2^(alpha) = 1.fraction
	unsigned short p = 256*doublesplittable.val;
	x = 256*doublesplittable.val/p - 1.0; // z := 256*2^(alpha)/p - 1
	double preroot = (1 + x*(1 + x*(x*(35 + x*(x*(4991 - 4055.1875*x)/8 - 805)/32)/8 - 7)/16)/8) * scalerootfactor[p-256]; // (1 + x)^(1/8) * (p/256)^(1/8)

	if( expnt < 0 ){ unsigned short expntdividedbyeight = (-expnt)/8;
/* Two Comments: -expnt <= 1022 = 127*8 + 6 ==> 0 <= expntdividedbyeight <= 127
 *               expnt = -8*expntdividedbyeight - (-8*expntdividedbyeight - expnt)
 *                    where (- expnt) - 8*expntdividedbyeight = (- expnt)%8;
 *                    then, 0 <=  (- expnt) - 8*expntdividedbyeight <= 7 */
		preroot /= eighthrootofpowersoftwo[(-expnt)-8*expntdividedbyeight];
		if(expntdividedbyeight < 64) return preroot/(1UL << expntdividedbyeight); // recall that the largest power of two fitting in an unsigned long is 2^(63)
		return (preroot/twoto64)/(1UL << (expntdividedbyeight - 64)); // expntdividedbyeight - 64 <= 127-64 = 63
	}
	if( expnt < 8 ) return preroot*eighthrootofpowersoftwo[expnt];
	return preroot*2.0*eighthrootofpowersoftwo[expnt-8]; // 8 <= expnt <= 15
}

/* BLAS-like convenience functions */
void VectorCopy_double(double *y, short n, double *x){ while(n--) *(y++) = *(x++); } /* Copies x in to y */
void VectorAddCntntTimesVector_double(double *y, short n, double a, double *x){ while(n--) *(y++) += a* *(x++); } /* Computes y = y + a * x */
void VectorSetToVectorPlusCntntTimesVector_double(double *y, short n, double *u, double a, double *x){ while(n--) *(y++) = *(u++) + a* *(x++); } /* Sets y = u + a * x */
void VectorSetToCntntTimesVector_double(double *y, short n, double a, double *x){ while(n--) *(y++) = a* *(x++); } /* Sets y = a * x */
void VectorMultiplyByCntnt_double(double *y, short n, double a){ while(n--) { *y = a* *y; y++; } } /* Sets y = a * y */
/* void VectorSetToCntntTimesVector_double(double *y, short n, double a, double *x){ while(n--) *(y++) = a* *(x++); } // Sets y = a * x
   void VectorAddVectorToMultiple_double(double *y, short n, double a, double *x){ while(n--) { *y = *(y) *a + *(x++); y++; } } // Computes y = x + a * y */
