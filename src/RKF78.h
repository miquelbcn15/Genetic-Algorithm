/* Runge-Kutta-Fehlberg 78 with adaptive stepsize version 1.1
 * implemented by Lluis Alseda on Jan 15, 2020
 * This version corrects the error:
 *     bij += 2; // bij = &b_{i+1,1} and (bij += 2) = &b_{i+1,3}
 * to
 *     bij += 3; // bij = &b_{i+1,1} and (bij += 3) = &b_{i+1,4} */
#ifndef _RKF78H_
#define _RKF78H_
#include <stdio.h>
#include <stdlib.h>
#define MIN(x,y) ((x) < (y) ? (x) : (y))

double RKF78(double *, double *,
             double *, double, double,
             double,
             void *,
             void (*)(double, double, double *, void *));

void RKF782tfin(double *, double *, double, void *,
             void (*)(double, double, double *, void *));
void ExitError(const char*, int);

void InitializeRKF78Sys(unsigned char);
double RKF78Sys(double *, double *,
             double *, double, double,
             double,
             void *,
             void (*)(double, double *, unsigned char, double *, void *));
void RKF78Sys2tfin(double *, double *, double, void *,
             void (*)(double, double *, unsigned char, double *, void *));
#endif
