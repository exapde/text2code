#include "forloop.h"
#include <cmath>

void forloop(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];

    dstype x2 = x0*x1;

    f[0 * N + i] = x2*uq0;
    f[1 * N + i] = x2*uq1;
    f[2 * N + i] = x2*uq2;
  }
}

void forloopjac(dstype* f, dstype* J1, dstype* J2, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];
    dstype x2 = x0*x1;

    f[0 * N + i] = x2*uq0;
    f[1 * N + i] = x2*uq1;
    f[2 * N + i] = x2*uq2;
    J1[0 * N + i] = x2;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = 0;
    J1[4 * N + i] = x2;
    J1[5 * N + i] = 0;
    J1[6 * N + i] = 0;
    J1[7 * N + i] = 0;
    J1[8 * N + i] = x2;
    J2[0 * N + i] = 0;
    J2[1 * N + i] = 0;
    J2[2 * N + i] = 0;
  }
}

void forloopjachess(dstype* f, dstype* J1, dstype* J2, dstype* H1, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];
    dstype x2 = x0*x1;

    f[0 * N + i] = x2*uq0;
    f[1 * N + i] = x2*uq1;
    f[2 * N + i] = x2*uq2;
    J1[0 * N + i] = x2;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = 0;
    J1[4 * N + i] = x2;
    J1[5 * N + i] = 0;
    J1[6 * N + i] = 0;
    J1[7 * N + i] = 0;
    J1[8 * N + i] = x2;
    J2[0 * N + i] = 0;
    J2[1 * N + i] = 0;
    J2[2 * N + i] = 0;
    H1[0 * N + i] = 0;
    H1[1 * N + i] = 0;
    H1[2 * N + i] = 0;
    H1[3 * N + i] = 0;
    H1[4 * N + i] = 0;
    H1[5 * N + i] = 0;
    H1[6 * N + i] = 0;
    H1[7 * N + i] = 0;
    H1[8 * N + i] = 0;
    H1[9 * N + i] = 0;
    H1[10 * N + i] = 0;
    H1[11 * N + i] = 0;
    H1[12 * N + i] = 0;
    H1[13 * N + i] = 0;
    H1[14 * N + i] = 0;
    H1[15 * N + i] = 0;
    H1[16 * N + i] = 0;
    H1[17 * N + i] = 0;
    H1[18 * N + i] = 0;
    H1[19 * N + i] = 0;
    H1[20 * N + i] = 0;
    H1[21 * N + i] = 0;
    H1[22 * N + i] = 0;
    H1[23 * N + i] = 0;
    H1[24 * N + i] = 0;
    H1[25 * N + i] = 0;
    H1[26 * N + i] = 0;
  }
}

