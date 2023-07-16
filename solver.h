#ifndef _Solver_
#define _Solver_

class Solver
{
private:
    // constants with default value
    const double GAMA = 1.4;
    const int KL = 2, KR = 3;  // ghost cells number of Left and Right
    int scheme = 3;
    int problem = 0;

public:
    double L = 1.0;      // domain length
    double TT = 0.1;           // Total time
    double CFL = 0.7;
    int J = 500;               // includes ghost cells
    double dt = 1e-6;
    double dx;
    double *d, *u, *p, *c;     // primitive
    double **fluxp, **fluxn, **flux; // flux
    double **Q, **Qn;          // conservative
    double *x;                 // location

    Solver();

    ~Solver();

    void init();

    void initShuOsher();

    void initSod();

    int computeDt();

    void computePrimitive();

    void splitFlux();

    void split_flux_LF_local();

    void split_flux_LF_global();

    void split_flux_SW();

    double diffP(double *fluxPositive);

    double diffN(double *fluxNegative);

    void computeFlux();

    void timeAdvance(int K);

    void q2qn();
};

#endif
