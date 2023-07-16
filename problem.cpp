#include "solver.h"
#include <cmath>

void Solver::init()
{
	d = new double[J];
	u = new double[J];
	p = new double[J];
	c = new double[J];
	x = new double[J];

	Q = new double *[J];
	Qn = new double *[J];
	fluxp = new double *[J];
	fluxn = new double *[J];
	flux = new double *[J];

	for (int i = 0; i < J; i++)
	{
		Q[i] = new double[3];
		Qn[i] = new double[3];
		fluxp[i] = new double[3];
		fluxn[i] = new double[3];
		flux[i] = new double[3];
	}

	for (int i = 0; i < J; i++)
		x[i] = i * dx - L / 2.0;

	switch (problem)
	{
	case 0:
		initSod();
		break;
	default:
		initShuOsher();
		break;
	}
}

void Solver::initShuOsher()
{
	double d1 = 3.857, u1 = 2.629, p1 = 10.333;
	double d2 = 0, u2 = 0.0, p2 = 1;

	for (int i = 0; i <= int(J * 1 / 10.0); i++)
	{
		Q[i][0] = d1;
		Q[i][1] = d1 * u1;
		Q[i][2] = p1 / (GAMA - 1) + d1 * u1 * u1 / 2;
	}
	for (int i = int(J * 1 / 10.0) + 1; i < J; i++)
	{
		d2 = 1 + 0.3 * sin(40 * x[i]);
		Q[i][0] = d2;
		Q[i][1] = d2 * u2;
		Q[i][2] = p2 / (GAMA - 1) + d2 * u2 * u2 / 2;
	}
}

void Solver::initSod()
{
	double d1 = 1.0, u1 = 0.0, p1 = 1.0;
	double d2 = 0.125, u2 = 0.0, p2 = 0.1;

	for (int i = 0; i <= int(J * 1 / 2.0); i++)
	{
		Q[i][0] = d1;
		Q[i][1] = d1 * u1;
		Q[i][2] = p1 / (GAMA - 1) + d1 * u1 * u1 / 2;
	}
	for (int i = int(J * 1 / 2.0) + 1; i < J; i++)
	{
		Q[i][0] = d2;
		Q[i][1] = d2 * u2;
		Q[i][2] = p2 / (GAMA - 1) + d2 * u2 * u2 / 2;
	}
}
