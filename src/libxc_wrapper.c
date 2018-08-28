/*
    Hong-Zhou Ye | Aug. 27, 2018

    A wrapper for libxc.

    The original purpose of this is for Tianyu's Fortran code on MPE, but I
    intentionally make it generic so might find broader use later.
*/


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "libxc_utils.h"


/*
    Evaluate xc energies and Exx given grid and MO coefficients.

    Input:
        np          number of grid points
        nbas        number of basis functions
        nocc        number of electron pairs
        nfunc       number of xc evaluations requested
        rdm1        [nbas*nbas] rdm1 in AO basis
        grid        [4*np] grid point positions (xyz) + weights
        ao_grid     [4*np*nbas] AO values and derivatives (xyz) on grid points
        func_ids    [nfunc] ID number of functionals
    Output:
        Es          [nfunc] nfunc elements correspond to func_ids
*/
int eval_xc_(const int * np_, const int * nbas_, 
	const int * nocc_, const int * nfunc_,
    const double * rdm1, const double * grid, const double * ao_grid,
    const int * func_ids, double * Es)
{
	// deference pointer
	const int np = *np_;
	const int nbas = *nbas_;
	const int nocc = *nocc_;
	const int nfunc = *nfunc_;

    const double * weight = grid + np*3;
    double * rho = (double *) malloc(np * sizeof(double));
    double * sigma = (double *) malloc(np * sizeof(double));

    // compute rho and sigma from grid and rdm1
    compute_rho(nbas, np, rdm1, ao_grid, rho);
    compute_sigma(nbas, np, rdm1, ao_grid, rho, sigma);

    // check that rho integrates to nelec = nocc*2
    if(check_nelec_grid(np, nocc, rho, weight) != 0)    exit(1);

    // compute xc energies
    double * exc = (double *) malloc(np * sizeof(double));
    for(int i = 0; i < nfunc; i++)
    {
        const int func_id = func_ids[i];
        if(compute_exc(np, func_id, rho, sigma, exc) != 0)  exit(1);
        Es[i] = integral_on_grid(np, rho, exc, weight);
    }

    free(rho);
    free(sigma);

    return 0;
}
