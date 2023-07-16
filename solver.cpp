#include "solver.h"
#include <iostream>
#include <cmath>

Solver::Solver()
{
	std ::cout << "Problem type: 0 for Sod, 1 for Shu-Osher" << std::endl;
	std::cin >> problem;
	std::cout << "Total time of simulation: " << std::endl;
	std::cin >> TT;
    std::cout << "Total length of domain: " << std::endl;
    std::cin >> L;
	std::cout << "Total bumber of points: " << std::endl;
	std::cin >> J;
	std::cout << "CFL number: " << std::endl;
	std::cin >> CFL;
	std ::cout << "Flux split schem: 1 for LF_local, 2 for Lf_global, 3 for SW" << std::endl;
	std::cin >> scheme;

	dx = L / (J - 1);
}

Solver::~Solver()
{
	delete[] d;
	delete[] u;
	delete[] p;
	delete[] c;
	delete[] x;
	for (int i = 0; i < J; i++)
	{
		delete[] Q[i];
		delete[] Qn[i];
		delete[] fluxp[i];
		delete[] fluxn[i];
		delete[] flux[i];
	}
	delete[] Q;
	delete[] Qn;
	delete[] fluxp;
	delete[] fluxn;
	delete[] flux;
	std::cout << "Deallocated variables, we were done" << std::endl;
}

int Solver::computeDt()
{
	double dtMin = 100;
	for (int i = 0; i < J; ++i)
	{
		double dtLocal = dx / (fabs(u[i]) + c[i]);
		if (dtLocal < dtMin)
			dtMin = dtLocal;
	}
	dt = CFL * dtMin;
	return 0;
}

void Solver::computePrimitive()
{
	for (int i = 0; i < J; i++)
	{
		d[i] = Q[i][0];
		u[i] = Q[i][1] / d[i];
		p[i] = (GAMA - 1.0) * (Q[i][2] - 0.5 * Q[i][1] * u[i]);
		c[i] = sqrt(GAMA * p[i] / d[i]);
	}
}

// 3rd Oder R-K
void Solver::timeAdvance(int KRK)
{
	double df = 0;
	for (int i = 1 + KL; i < J - KR; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			df = -(flux[i][j] - flux[i - 1][j]) / dx;
			if (KRK == 1)
			{
				Q[i][j] = Qn[i][j] + dt * df;
			}
			else if (KRK == 2)
			{
				Q[i][j] = 0.75 * Qn[i][j] + 0.25 * (Q[i][j] + dt * df);
			}
			else if (KRK == 3)
			{
				Q[i][j] = 1.0 / 3.0 * Qn[i][j] + 2.0 / 3.0 * (Q[i][j] + dt * df);
			}
		}
	}
}

void Solver::q2qn()
{
	for (int i = 0; i < J; i++)
	{
		for (int j = 0; j < 3; j++)
			Qn[i][j] = Q[i][j];
	}
}
