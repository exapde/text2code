#include "mass.h"
#include <cmath>

void mass(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];


    f[0 * N + i] = pow(uq2, 2)*pow(uq1, 2) - pow(uq0, 2);
  }
}

void massjac(dstype* f, dstype* J1, dstype* J2, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype x0 = pow(uq2, 2);
    dstype x1 = pow(uq1, 2);

    f[0 * N + i] = x0*x1 - pow(uq0, 2);
    J1[0 * N + i] = -2*uq0;
    J1[1 * N + i] = 2*x0*uq1;
    J1[2 * N + i] = 2*x1*uq2;
    J2[0 * N + i] = 0;
  }
}

void massjachess(dstype* f, dstype* J1, dstype* J2, dstype* H1, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype x0 = pow(uq2, 2);
    dstype x1 = pow(uq1, 2);
    dstype x2 = 2*x0;
    dstype x3 = 2*x1;
    dstype x4 = 4*uq2*uq1;

    f[0 * N + i] = x0*x1 - pow(uq0, 2);
    J1[0 * N + i] = -2*uq0;
    J1[1 * N + i] = x2*uq1;
    J1[2 * N + i] = x3*uq2;
    J2[0 * N + i] = 0;
    H1[0 * N + i] = -2;
    H1[1 * N + i] = 0;
    H1[2 * N + i] = 0;
    H1[3 * N + i] = 0;
    H1[4 * N + i] = x2;
    H1[5 * N + i] = x4;
    H1[6 * N + i] = 0;
    H1[7 * N + i] = x4;
    H1[8 * N + i] = x3;
  }
}

