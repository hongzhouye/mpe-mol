#ifndef _LIBXC_UTILS_H_
#define _LIBXC_UTILS_H_

#include <xc.h>

#define LDA_X 1
#define LDA_C_VWN 7
#define GGA_X_PBE 101
#define GGA_C_PBE 130
#define GGA_X_B88 106
#define GGA_C_LYP 131

void compute_rho(int nbas, int np, const double * rdm1, const double * ao_grid,
    double * rho);
void compute_sigma(int nbas, int np, const double * rdm1,
    const double * ao_grid, const double * rho, double * sigma);
int check_nelec_grid(int np, int nocc, const double * rho,
    const double * weight);
double integral_on_grid(int np, const double * rho, const double * ex,
    const double * weight);
int compute_exc(const int np, const int func_id,
    const double * rho, const double * sigma, double * exc);
double compute_EJ(int nbas, int nocc, const double * Vmo);
double compute_exx(int nbas, int nocc, const double * Vmo);

#endif
