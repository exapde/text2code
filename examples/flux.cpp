#include "flux.h"
#include <cmath>

void flux(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype mu0 = mu[0*N+i];

    dstype x0 = mu0*(2 + pow(uq0, 2) + tanh(uq0));

    f[0 * N + i] = x0*uq1;
    f[1 * N + i] = x0*uq2;
  }
}

void fluxjac(dstype* f, dstype* J1, dstype* J2, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype mu0 = mu[0*N+i];
    dstype x0 = tanh(uq0);
    dstype x1 = mu0*(2 + x0 + pow(uq0, 2));
    dstype x2 = mu0*(1 + 2*uq0 - pow(x0, 2));

    f[0 * N + i] = x1*uq1;
    f[1 * N + i] = x1*uq2;
    J1[0 * N + i] = x2*uq1;
    J1[1 * N + i] = x1;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = x2*uq2;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = x1;
    J2[0 * N + i] = 0;
    J2[1 * N + i] = 0;
  }
}

void fluxjachess(dstype* f, dstype* J1, dstype* J2, dstype* H1, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype mu0 = mu[0*N+i];
    dstype x0 = tanh(uq0);
    dstype x1 = mu0*(2 + x0 + pow(uq0, 2));
    dstype x2 = 1 - pow(x0, 2);
    dstype x3 = mu0*(2*uq0 + x2);
    dstype x4 = mu0*(2 - 2*x0*x2);

    f[0 * N + i] = x1*uq1;
    f[1 * N + i] = x1*uq2;
    J1[0 * N + i] = x3*uq1;
    J1[1 * N + i] = x1;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = x3*uq2;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = x1;
    J2[0 * N + i] = 0;
    J2[1 * N + i] = 0;
    H1[0 * N + i] = x4*uq1;
    H1[1 * N + i] = x3;
    H1[2 * N + i] = 0;
    H1[3 * N + i] = x3;
    H1[4 * N + i] = 0;
    H1[5 * N + i] = 0;
    H1[6 * N + i] = 0;
    H1[7 * N + i] = 0;
    H1[8 * N + i] = 0;
    H1[9 * N + i] = x4*uq2;
    H1[10 * N + i] = 0;
    H1[11 * N + i] = x3;
    H1[12 * N + i] = 0;
    H1[13 * N + i] = 0;
    H1[14 * N + i] = 0;
    H1[15 * N + i] = x3;
    H1[16 * N + i] = 0;
    H1[17 * N + i] = 0;
  }
}

