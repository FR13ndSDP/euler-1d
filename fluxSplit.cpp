#include "solver.h"
#include "cmath"

void Solver::splitFlux()
{
	switch (scheme)
	{
	case 0:
		return split_flux_LF_local();
		break;
	case 1:
		return split_flux_LF_global();
		break;
	default:
		return split_flux_SW();
		break;
	}
}

// local L-F flux splitting
void Solver::split_flux_LF_local()
{
	// local max eigenvalue
	double lmax;
	for (int i = 0; i < J; i++)
	{
		// lmax = fmax(fmax(fabs(u[i]), fabs(u[i] - c[i])), fabs(u[i] + c[i]));
		lmax = fabs(u[i]) + c[i];
		fluxp[i][0] = 0.50 * (d[i] * u[i] + lmax * Q[i][0]);
		fluxp[i][1] = 0.50 * (d[i] * u[i] * u[i] + p[i] + lmax * Q[i][1]);
		fluxp[i][2] = 0.50 * ((Q[i][2] + p[i]) * u[i] + lmax * Q[i][2]);
		fluxn[i][0] = 0.50 * (d[i] * u[i] - lmax * Q[i][0]);
		fluxn[i][1] = 0.50 * (d[i] * u[i] * u[i] + p[i] - lmax * Q[i][1]);
		fluxn[i][2] = 0.50 * ((Q[i][2] + p[i]) * u[i] - lmax * Q[i][2]);
	}
}

void Solver::split_flux_LF_global()
{
	double lmax, gmax = 0;
	for (int i = 0; i < J; i++)
	{
		lmax = fabs(u[i]) + c[i];
		// global max eigenvalue
		gmax = fmax(gmax, lmax);
	}
	for (int i = 0; i < J; i++)
	{
		fluxp[i][0] = 0.50 * (d[i] * u[i] + gmax * Q[i][0]);
		fluxp[i][1] = 0.50 * (d[i] * u[i] * u[i] + p[i] + gmax * Q[i][1]);
		fluxp[i][2] = 0.50 * ((Q[i][2] + p[i]) * u[i] + gmax * Q[i][2]);
		fluxn[i][0] = 0.50 * (d[i] * u[i] - gmax * Q[i][0]);
		fluxn[i][1] = 0.50 * (d[i] * u[i] * u[i] + p[i] - gmax * Q[i][1]);
		fluxn[i][2] = 0.50 * ((Q[i][2] + p[i]) * u[i] - gmax * Q[i][2]);
	}
}

void Solver::split_flux_SW()
{
	double El1, El2, El3, El1p, El1m, El2p, El2m, El3p, El3m, tmp, WP, WM;
	double constVar1 = GAMA - 1.0;
	double constVar2 = 2.0 * constVar1;
	double constVar3 = (3.0 - GAMA) / constVar2;
	for (int i = 0; i < J; i++)
	{
		El1 = u[i];
		El2 = u[i] - c[i];
		El3 = u[i] + c[i];
		El1p = 0.50 * (El1 + fabs(El1));
		El1m = 0.50 * (El1 - fabs(El1));
		El2p = 0.50 * (El2 + fabs(El2));
		El2m = 0.50 * (El2 - fabs(El2));
		El3p = 0.50 * (El3 + fabs(El3));
		El3m = 0.50 * (El3 - fabs(El3));
		tmp = d[i] / (2.0 * GAMA);

		fluxp[i][0] = tmp * (constVar2 * El1p + El2p + El3p);
		fluxn[i][0] = tmp * (constVar2 * El1m + El2m + El3m);
		fluxp[i][1] = tmp * (constVar2 * El1p * El1 + El2p * El2 + El3p * El3);
		fluxn[i][1] = tmp * (constVar2 * El1m * El1 + El2m * El2 + El3m * El3);
		WP = (El2p + El3p) * c[i] * c[i] * constVar3;
		WM = (El2m + El3m) * c[i] * c[i] * constVar3;
		fluxp[i][2] = tmp * (constVar1 * El1p * u[i] * u[i]\
		              + El2p / 2.0 * pow(El2, 2)\
		              + WP\
					  + El3p / 2.0 * pow(El3, 2));
		fluxn[i][2] = tmp * (constVar1 * El1m * u[i] * u[i]\
		              + El2m / 2.0 * pow(El2, 2)\
		              + WM\
					  + El3m / 2.0 * pow(El3, 2));
	}
}