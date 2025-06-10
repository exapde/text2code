#include "SymbolicFunctions.hpp"

std::vector<Expression> conductivity(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta) {
    std::vector<Expression> kappa;
    kappa.resize(1);

    kappa[0]  =  mu[0]*(2 + uq[0]*uq[0] + tanh(uq[0]));
    return kappa;
}

std::vector<Expression> flux(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta) {
    std::vector<Expression> f;
    f.resize(2);

    auto kappa = conductivity(uq, w, v, x, time, mu, eta);
    f[0]  =  kappa[0]*uq[1];
    f[1]  =  kappa[0]*uq[2];
    return f;
}

std::vector<Expression> forloop(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta) {
    std::vector<Expression> f;
    f.resize(3);

    for (int i = 0; i <= 2; ++i) {
    f[i]  =  uq[i] * x[0] * x[1];
    }
    return f;
}

std::vector<Expression> mass(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta) {
    std::vector<Expression> f;
    f.resize(1);

    SymEngine::DenseMatrix A(2, 2);
    A.set(0, 0, uq[1]*uq[1]);
    A.set(0, 1, uq[0]);
    A.set(1, 0, uq[0]);
    A.set(1, 1, uq[2]*uq[2]);
    f[0] = A.det();
    return f;
}

std::vector<Expression> source(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta) {
    std::vector<Expression> f;
    f.resize(1);

    Expression x1 = x[0];
    Expression x2 = x[1];
    f[0]  =  x1*sin(x2)*sin(SymEngine::pi*time);
    return f;
}

std::vector<Expression> fbouhdg(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau) {
    std::vector<Expression> fb;
    fb.resize(2);

    fb[0]  =  tau[0]*(mu[1] - uhat[0]);
    auto f = flux(uq, w, v, x, time, mu, eta);
    fb[1]  =  f[0]*n[0] + f[1]*n[1] + tau[0]*(uq[0]-uhat[0]);
    return fb;
}

std::vector<Expression> initu(const std::vector<Expression>& x) {
    std::vector<Expression> f;
    f.resize(2);

    f[0]  =  x[0]*x[0] + x[1];
    f[1]  =  x[0]*x[0] - x[1];
    return f;
}

