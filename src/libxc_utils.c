#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "libxc_utils.h"

const double thresh = 1.E-2;

const int known_func_ids[6] = {
    1,      // LDA_X
    7,      // LDA_C_VWN
    101,    // GGA_X_PBE
    130,    // GGA_C_PBE
    106,    // GGA_X_B88
    131     // GGA_C_LYP
};

void compute_rho(int nbas, int np, const double * rdm1, const double * ao_grid,
    double * rho)
{
    for(int p = 0; p < np; p++)
    {
        double temp = 0.;
        for(int mu = 0; mu < nbas; mu++)
        {
            const double * phi_mu = ao_grid + mu*np;
            for(int nu = 0; nu < nbas; nu++)
            {
                const double * phi_nu = ao_grid + nu*np;
                temp += rdm1[mu*nbas+nu] * phi_mu[p] * phi_nu[p];
            }
        }
        rho[p] = temp;
    }
}

void compute_sigma(int nbas, int np, const double * rdm1,
    const double * ao_grid, const double * rho, double * sigma)
{
    const double exponent = 4. / 3.;
    const double factor = 0.5 / pow(3.*M_PI*M_PI, 0.25*exponent);
    for(int p = 0; p < np; p++)
    {
        double temp = 0;
        for(int alp = 1; alp < 4; alp++)
        {
            double temp_alp = 0.;
            const double * dphidalp = ao_grid + alp*nbas*np;
            for(int mu = 0; mu < nbas; mu++)
            {
                const double * phi_mu = ao_grid + mu*np;
                for(int nu = 0; nu < nbas; nu++)
                {
                    const double * dphidalp_nu = dphidalp + nu*np;
                    temp_alp += rdm1[mu*nbas+nu] * phi_mu[p] * dphidalp_nu[p];
                }
            }
            temp_alp *= 2.;
            temp += temp_alp * temp_alp;
        }
        // sigma[p] = factor * sqrt(temp) / pow(rho[p], exponent);
        // sigma[p] = sqrt(temp);
        sigma[p] = temp;
    }
}

int check_nelec_grid(int np, int nocc, const double * rho,
    const double * weight)
{
    double nelec = 0.;
    for(int p = 0; p < np; p++)
        nelec += rho[p] * weight[p];

    if(fabs(nelec-nocc*2.) > thresh)
    {
        fprintf(stderr, "Computed total number of electrons (%.5f) is not correct\n", nelec);
        return 1;
    }

    return 0;
}

double integral_on_grid(int np, const double * rho, const double * ex,
    const double * weight)
{
    double sum = 0.;
    for(int p = 0; p < np; p++)
        sum += rho[p] * ex[p] * weight[p];

    return sum;
}

int compute_exc(const int np, const int func_id,
    const double * rho, const double * sigma, double * exc)
{
    // check whether the input func_id is known
    bool known_id = false;
    for(int i = 0; i < sizeof(known_func_ids)/sizeof(known_func_ids[0]); i++)
        if(func_id == known_func_ids[i])
        {
            known_id = true;
            break;
        }
    if(!known_id)
    {
        fprintf(stderr, "Input func '%d' is not supported.\n", func_id);
        return 1;
    }

    // initialize libxc func object
    xc_func_type func;
    if(xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0)
    {
        fprintf(stderr, "Failed to initialize xc func '%d'.\n", func_id);
        return 1;
    }

    switch(func.info->family)
    {
        case XC_FAMILY_LDA:
            xc_lda_exc(&func, np, rho, exc);
            break;
        case XC_FAMILY_GGA:
        case XC_FAMILY_HYB_GGA:
            xc_gga_exc(&func, np, rho, sigma, exc);
            break;
    }

    printf("Computing %s energy\n", xc_functional_get_name(func_id));

    xc_func_end(&func);

    return 0;
}

#define IND4(i,j,k,l,n) i*n*n*n+j*n*n+k*n+l
double compute_EJ(int nbas, int nocc, const double * Vmo)
{
    double EJ = 0.;
    for(int i = 0; i < nocc; i++)   for(int j = 0; j < nocc; j++)
        EJ += Vmo[IND4(i,i,j,j,nbas)];
    EJ *= 2.;

    return EJ;
}

double compute_exx(int nbas, int nocc, const double * Vmo)
{
    double Exx = 0.;
    for(int i = 0; i < nocc; i++)   for(int j = 0; j < nocc; j++)
        Exx -= Vmo[IND4(i,j,j,i,nbas)];

    return Exx;
}
