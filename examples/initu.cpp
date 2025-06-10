#include "initu.h"
#include <cmath>

void initu(dstype* f, const dstype* x, const int N, const int szf, const int szx)
{

  for (int i = 0; i < N; ++i) {
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];

    dstype x2 = pow(x0, 2);

    f[0 * N + i] = x1 + x2;
    f[1 * N + i] = -x1 + x2;
  }
}

