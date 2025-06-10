#include "fbouhdg.h"
#include <cmath>

void fbouhdg(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const dstype* uhat, const dstype* n, const dstype* tau, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta, const int szuhat, const int szn, const int sztau)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype mu0 = mu[0*N+i];
    dstype mu1 = mu[1*N+i];
    dstype uhat0 = uhat[0*N+i];
    dstype n0 = n[0*N+i];
    dstype n1 = n[1*N+i];
    dstype tau0 = tau[0*N+i];

    dstype x0 = -uhat0;
    dstype x1 = mu0*(2 + pow(uq0, 2) + tanh(uq0));

    f[0 * N + i] = tau0*(mu1 + x0);
    f[1 * N + i] = tau0*(uq0 + x0) + n0*x1*uq1 + n1*x1*uq2;
  }
}

void fbouhdgjac(dstype* f, dstype* J1, dstype* J2, dstype* J3, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const dstype* uhat, const dstype* n, const dstype* tau, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta, const int szuhat, const int szn, const int sztau)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype mu0 = mu[0*N+i];
    dstype mu1 = mu[1*N+i];
    dstype uhat0 = uhat[0*N+i];
    dstype n0 = n[0*N+i];
    dstype n1 = n[1*N+i];
    dstype tau0 = tau[0*N+i];
    dstype x0 = -uhat0;
    dstype x1 = tanh(uq0);
    dstype x2 = mu0*(2 + x1 + pow(uq0, 2));
    dstype x3 = n0*x2;
    dstype x4 = n1*x2;
    dstype x5 = mu0*(1 + 2*uq0 - pow(x1, 2));
    dstype x6 = -tau0;

    f[0 * N + i] = tau0*(mu1 + x0);
    f[1 * N + i] = tau0*(uq0 + x0) + x3*uq1 + x4*uq2;
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = tau0 + n0*x5*uq1 + n1*x5*uq2;
    J1[4 * N + i] = x3;
    J1[5 * N + i] = x4;
    J2[0 * N + i] = 0;
    J2[1 * N + i] = 0;
    J3[0 * N + i] = x6;
    J3[1 * N + i] = x6;
  }
}

void fbouhdgjachess(dstype* f, dstype* J1, dstype* J2, dstype* J3, dstype* H1, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const dstype* uhat, const dstype* n, const dstype* tau, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta, const int szuhat, const int szn, const int sztau)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype mu0 = mu[0*N+i];
    dstype mu1 = mu[1*N+i];
    dstype uhat0 = uhat[0*N+i];
    dstype n0 = n[0*N+i];
    dstype n1 = n[1*N+i];
    dstype tau0 = tau[0*N+i];
    dstype x0 = -uhat0;
    dstype x1 = tanh(uq0);
    dstype x2 = 2 + x1 + pow(uq0, 2);
    dstype x3 = n0*mu0;
    dstype x4 = x2*x3;
    dstype x5 = n1*mu0;
    dstype x6 = x2*x5;
    dstype x7 = 1 - pow(x1, 2);
    dstype x8 = 2*uq0 + x7;
    dstype x9 = x5*x8;
    dstype x10 = x3*x8;
    dstype x11 = -tau0;
    dstype x12 = 2 - 2*x1*x7;

    f[0 * N + i] = tau0*(mu1 + x0);
    f[1 * N + i] = tau0*(uq0 + x0) + x4*uq1 + x6*uq2;
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = tau0 + uq1*x10 + x9*uq2;
    J1[4 * N + i] = x4;
    J1[5 * N + i] = x6;
    J2[0 * N + i] = 0;
    J2[1 * N + i] = 0;
    J3[0 * N + i] = x11;
    J3[1 * N + i] = x11;
    H1[0 * N + i] = 0;
    H1[1 * N + i] = 0;
    H1[2 * N + i] = 0;
    H1[3 * N + i] = 0;
    H1[4 * N + i] = 0;
    H1[5 * N + i] = 0;
    H1[6 * N + i] = 0;
    H1[7 * N + i] = 0;
    H1[8 * N + i] = 0;
    H1[9 * N + i] = x3*uq1*x12 + x5*uq2*x12;
    H1[10 * N + i] = x10;
    H1[11 * N + i] = x9;
    H1[12 * N + i] = x10;
    H1[13 * N + i] = 0;
    H1[14 * N + i] = 0;
    H1[15 * N + i] = x9;
    H1[16 * N + i] = 0;
    H1[17 * N + i] = 0;
  }
}

