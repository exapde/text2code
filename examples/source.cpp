#include "source.h"
#include <cmath>

void source(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];


    f[0 * N + i] = x0*sin(x1)*sin(acos(-1)*time);
  }
}

void sourcejac(dstype* f, dstype* J1, dstype* J2, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];

    f[0 * N + i] = x0*sin(x1)*sin(acos(-1)*time);
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J2[0 * N + i] = 0;
  }
}

void sourcejachess(dstype* f, dstype* J1, dstype* J2, dstype* H1, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)
{

  for (int i = 0; i < N; ++i) {
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];

    f[0 * N + i] = x0*sin(x1)*sin(acos(-1)*time);
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J2[0 * N + i] = 0;
    H1[0 * N + i] = 0;
    H1[1 * N + i] = 0;
    H1[2 * N + i] = 0;
    H1[3 * N + i] = 0;
    H1[4 * N + i] = 0;
    H1[5 * N + i] = 0;
    H1[6 * N + i] = 0;
    H1[7 * N + i] = 0;
    H1[8 * N + i] = 0;
  }
}

