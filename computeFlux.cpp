#include "solver.h"
#include <cmath>

double Solver::diffP(double *fluxLocal)
{
	double fpx;
	double eps = 1e-6, p = 2;
	double c1 = 1.0 / 10.0, c2 = 3.0 / 5.0, c3 = 3.0 / 10.0;

	double f1 = 1.0 / 3.0 * fluxLocal[KL - 2] \
	            - 7.0 / 6.0 * fluxLocal[KL - 1]\
	            + 11.0 / 6.0 * fluxLocal[KL];
	double f2 = -1.0 / 6.0 * fluxLocal[KL - 1]\
	            + 5.0 / 6.0 * fluxLocal[KL] \
	            + 1.0 / 3.0 * fluxLocal[KL+1];
	double f3 = 1.0 / 3.0 * fluxLocal[KL]\
	            + 5.0 / 6.0 * fluxLocal[KL + 1]\
				- 1.0 / 6.0 * fluxLocal[KL + 2];

	double IS1 = 1.0 / 4.0 * pow(fluxLocal[KL - 2] - 4 * fluxLocal[KL - 1] + 3 * fluxLocal[KL], 2) \
	             + 13.0 / 12.0 * pow(fluxLocal[KL - 2] - 2 * fluxLocal[KL - 1] + fluxLocal[KL], 2);
	double IS2 = 1.0 / 4.0 * pow(fluxLocal[KL - 1] - fluxLocal[KL + 1], 2)\
	             + 13.0 / 12.0 * pow(fluxLocal[KL - 1] - 2 * fluxLocal[KL] + fluxLocal[KL + 1], 2);
	double IS3 = 1.0 / 4.0 * pow(3 * fluxLocal[KL] - 4 * fluxLocal[KL + 1] + fluxLocal[KL + 2], 2)\
	             + 13.0 / 12.0 * pow(fluxLocal[KL] - 2 * fluxLocal[KL + 1] + fluxLocal[KL + 2], 2);

	double alpha1 = c1 / pow(eps + IS1, p);
	double alpha2 = c2 / pow(eps + IS2, p);
	double alpha3 = c3 / pow(eps + IS3, p);
	double omega1 = alpha1 / (alpha1 + alpha2 + alpha3);
	double omega2 = alpha2 / (alpha1 + alpha2 + alpha3);
	double omega3 = alpha3 / (alpha1 + alpha2 + alpha3);

	fpx = omega1 * f1 + omega2 * f2 + omega3 * f3;
	return fpx;
}

double Solver::diffN(double *fluxLocal)
{
	double fnx;
	double eps = 1e-6, p = 2;
	double c1 = 1.0 / 10.0, c2 = 3.0 / 5.0, c3 = 3.0 / 10.0;

	double f1 = 1.0 / 3.0 * fluxLocal[KL + 3] \
	            - 7.0 / 6.0 * fluxLocal[KL + 2]\
	            + 11.0 / 6.0 * fluxLocal[KL + 1];
	double f2 = -1.0 / 6.0 * fluxLocal[KL + 2]\
	            + 5.0 / 6.0 * fluxLocal[KL+ 1] \
	            + 1.0 / 3.0 * fluxLocal[KL];
	double f3 = 1.0 / 3.0 * fluxLocal[KL + 1]\
	            + 5.0 / 6.0 * fluxLocal[KL]\
				- 1.0 / 6.0 * fluxLocal[KL - 1];

	double IS1 = 1.0 / 4.0 * pow(fluxLocal[KL+3] - 4 * fluxLocal[KL+2] + 3 * fluxLocal[KL+1], 2) \
	             + 13.0 / 12.0 * pow(fluxLocal[KL+3] - 2 * fluxLocal[KL+2] + fluxLocal[KL + 1], 2);
	double IS2 = 1.0 / 4.0 * pow(fluxLocal[KL+2] - fluxLocal[KL], 2)\
	             + 13.0 / 12.0 * pow(fluxLocal[KL+2] - 2 * fluxLocal[KL+1] + fluxLocal[KL], 2);
	double IS3 = 1.0 / 4.0 * pow(3 * fluxLocal[KL+1] - 4 * fluxLocal[KL] + fluxLocal[KL-1], 2)\
	             + 13.0 / 12.0 * pow(fluxLocal[KL+1] - 2 * fluxLocal[KL] + fluxLocal[KL-1], 2);

	double alpha1 = c1 / pow(eps + IS1, p);
	double alpha2 = c2 / pow(eps + IS2, p);
	double alpha3 = c3 / pow(eps + IS3, p);
	double omega1 = alpha1 / (alpha1 + alpha2 + alpha3);
	double omega2 = alpha2 / (alpha1 + alpha2 + alpha3);
	double omega3 = alpha3 / (alpha1 + alpha2 + alpha3);

	fnx = omega1 * f1 + omega2 * f2 + omega3 * f3;
	return fnx;
}

void Solver::computeFlux()
{
	double *fluxLocalp = new double[KL + KR + 1], *fluxLocaln = new double[KL + KR + 1];
	double fpx = 0, fnx = 0;

	for (int i = KL; i < J - KR; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = -KL; k <= KR; k++)
			{
				{
					fluxLocalp[k + KL] = fluxp[i + k][j];
					fluxLocaln[k + KL] = fluxn[i + k][j];
				}
			}

			// 5th order WENO
			fpx = diffP(fluxLocalp);
			fnx = diffN(fluxLocaln);
			// flux in j+1/2
			flux[i][j] = fpx + fnx;
		}
	}

	delete[] fluxLocalp;
	delete[] fluxLocaln;
}