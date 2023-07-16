#include <cstdio>
#include <iostream>
#include "solver.h"
#ifdef TIMING
#include <time.h>
#endif

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cout << "No file name provided, exit" << std::endl;
		exit(1);
	}

	Solver euler;
	double T = 0;
	euler.init();
#ifdef TIMING
	clock_t start_t, end_t;
	start_t = clock();
#endif
	while (T < euler.TT)
	{
		euler.q2qn();
        for (auto KRK : {1, 2, 3})
		{
			euler.computePrimitive();
			euler.splitFlux();
			euler.computeFlux();
			euler.timeAdvance(KRK);
			if (KRK == 3)
			{
				T += euler.dt;
				printf("T=%10g\t dt=%10g\n", T, euler.dt);
			}
		}
		euler.computeDt();
		if (euler.dt + T > euler.TT)
			euler.dt = euler.TT - T;
	}
	euler.computePrimitive();
	FILE *fp;
	fp = fopen(argv[1], "w");
	for (int i = 0; i < euler.J; i++)
	{
		fprintf(fp, "%f %f %f %f\n", euler.x[i], euler.d[i], euler.u[i], euler.p[i]);
	}
	fclose(fp);
#ifdef TIMING
	end_t = clock();
	std ::cout << double(end_t-start_t)/CLOCKS_PER_SEC << " s elapsed" << std::endl;
#endif
	return 0;
}
