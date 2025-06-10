#include "SymbolicScalarsVectors.hpp"

SymbolicScalarsVectors::SymbolicScalarsVectors() {

    time =   Expression("time");
    a =   Expression("a");
    b =   Expression("b");
    c =   Expression("c");
    d =   Expression("d");
    e =   Expression("e");

    sztau = 1;
    tau.resize(1);
    for (int i = 0; i < 1; ++i) {
         tau[i] = Expression("tau"  + std::to_string(i));
    }

    szw = 1;
    w.resize(1);
    for (int i = 0; i < 1; ++i) {
         w[i] = Expression("w"  + std::to_string(i));
    }

    szeta = 0;
    eta.resize(0);
    for (int i = 0; i < 0; ++i) {
         eta[i] = Expression("eta"  + std::to_string(i));
    }

    szuq = 3;
    uq.resize(3);
    for (int i = 0; i < 3; ++i) {
         uq[i] = Expression("uq"  + std::to_string(i));
    }

    szn = 2;
    n.resize(2);
    for (int i = 0; i < 2; ++i) {
         n[i] = Expression("n"  + std::to_string(i));
    }

    szuhat = 1;
    uhat.resize(1);
    for (int i = 0; i < 1; ++i) {
         uhat[i] = Expression("uhat"  + std::to_string(i));
    }

    szv = 0;
    v.resize(0);
    for (int i = 0; i < 0; ++i) {
         v[i] = Expression("v"  + std::to_string(i));
    }

    szmu = 2;
    mu.resize(2);
    for (int i = 0; i < 2; ++i) {
         mu[i] = Expression("mu"  + std::to_string(i));
    }

    szx = 2;
    x.resize(2);
    for (int i = 0; i < 2; ++i) {
         x[i] = Expression("x"  + std::to_string(i));
    }

    outputfunctions.assign(7, false);
    outputfunctions[1] = true;
    outputfunctions[2] = true;
    outputfunctions[3] = true;
    outputfunctions[4] = true;
    outputfunctions[5] = true;
    outputfunctions[6] = true;

    funcnames = {"conductivity", "flux", "forloop", "mass", "source", "fbouhdg", "initu"};

    funcargs = {
        {"uq", "w", "v", "x", "time", "mu", "eta"},
        {"uq", "w", "v", "x", "time", "mu", "eta"},
        {"uq", "w", "v", "x", "time", "mu", "eta"},
        {"uq", "w", "v", "x", "time", "mu", "eta"},
        {"uq", "w", "v", "x", "time", "mu", "eta"},
        {"uq", "w", "v", "x", "time", "mu", "eta", "uhat", "n", "tau"},
        {"x"}
    };

    funcargssizes = {
        {"szuq", "szw", "szv", "szx", "sztime", "szmu", "szeta"},
        {"szuq", "szw", "szv", "szx", "sztime", "szmu", "szeta"},
        {"szuq", "szw", "szv", "szx", "sztime", "szmu", "szeta"},
        {"szuq", "szw", "szv", "szx", "sztime", "szmu", "szeta"},
        {"szuq", "szw", "szv", "szx", "sztime", "szmu", "szeta"},
        {"szuq", "szw", "szv", "szx", "sztime", "szmu", "szeta", "szuhat", "szn", "sztau"},
        {"szx"}
    };

    funcdecls = {
       "void conductivity(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "void flux(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "void forloop(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "void mass(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "void source(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "void fbouhdg(dstype* f, const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const dstype* uhat, const dstype* n, const dstype* tau, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta, const int szuhat, const int szn, const int sztau)", 
       "void initu(dstype* f, const dstype* x, const int N, const int szf, const int szx)"
    };

    funcjacdecls = {
       "const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta)", 
       "const dstype* uq, const dstype* w, const dstype* v, const dstype* x, const dstype time, const dstype* mu, const dstype* eta, const dstype* uhat, const dstype* n, const dstype* tau, const int N, const int szf, const int szuq, const int szw, const int szv, const int szx, const int szmu, const int szeta, const int szuhat, const int szn, const int sztau)", 
       "const dstype* x, const int N, const int szf, const int szx)"
    };

    inputvectors = {
        {{"uq", uq}, {"w", w}, {"v", v}, {"x", x}, {"mu", mu}, {"eta", eta}},
        {{"uq", uq}, {"w", w}, {"v", v}, {"x", x}, {"mu", mu}, {"eta", eta}},
        {{"uq", uq}, {"w", w}, {"v", v}, {"x", x}, {"mu", mu}, {"eta", eta}},
        {{"uq", uq}, {"w", w}, {"v", v}, {"x", x}, {"mu", mu}, {"eta", eta}},
        {{"uq", uq}, {"w", w}, {"v", v}, {"x", x}, {"mu", mu}, {"eta", eta}},
        {{"uq", uq}, {"w", w}, {"v", v}, {"x", x}, {"mu", mu}, {"eta", eta}, {"uhat", uhat}, {"n", n}, {"tau", tau}},
        {{"x", x}}
    };

    inputscalars = {
        {{"time", time}},
        {{"time", time}},
        {{"time", time}},
        {{"time", time}},
        {{"time", time}},
        {{"time", time}},
        {}
    };

    jacobianInputs = {
        {uq, w},
        {uq, w},
        {uq, w},
        {uq, w},
        {uq, w},
        {uq, w, uhat},
        {}
    };

    hessianInputs = {
        {uq},
        {uq},
        {uq},
        {uq},
        {uq},
        {uq},
        {}
    };

}

std::vector<Expression> SymbolicScalarsVectors::evaluateSymbolicFunctions(int call)
{
  std::vector<Expression> f;

  switch (call) {
    case 0:
      f = conductivity(uq, w, v, x, time, mu, eta);
      break;
    case 1:
      f = flux(uq, w, v, x, time, mu, eta);
      break;
    case 2:
      f = forloop(uq, w, v, x, time, mu, eta);
      break;
    case 3:
      f = mass(uq, w, v, x, time, mu, eta);
      break;
    case 4:
      f = source(uq, w, v, x, time, mu, eta);
      break;
    case 5:
      f = fbouhdg(uq, w, v, x, time, mu, eta, uhat, n, tau);
      break;
    case 6:
      f = initu(x);
      break;
    default:
      throw std::runtime_error("Invalid function call in evaluateSymbolicFunctions");
  }

  return f;
}

void SymbolicScalarsVectors::func2cse(vec_pair &replacements, vec_basic &reduced_exprs, const std::vector<Expression> &f) {

   vec_basic exprs;
   for (const auto &fi : f) {
       exprs.push_back(fi.get_basic());
   }
   cse(replacements, reduced_exprs, exprs);
}

void SymbolicScalarsVectors::funcjac2cse(vec_pair &replacements,
                                         vec_basic &reduced_exprs_f,
                                         std::vector<vec_basic> &reduced_exprs_J,
                                         const std::vector<Expression> &f,
                                         const std::vector<std::vector<Expression>>& inputs_J) {
    vec_basic exprs;

    // Track original sizes
    int n_f = static_cast<int>(f.size());
    std::vector<int> jacobian_sizes;

    // Add original function expressions
    for (const auto &fi : f) {
        exprs.push_back(fi.get_basic());
    }

    // Append Jacobians and record sizes for each input group
    for (const auto& input_vec : inputs_J) {
        int m = f.size();
        int n = input_vec.size();
        jacobian_sizes.push_back(m * n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                exprs.push_back(f[i].diff(input_vec[j]).get_basic());
            }
        }
    }

    // Apply CSE
    vec_basic reduced_exprs;
    cse(replacements, reduced_exprs, exprs);

    // Decompose: f
    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + n_f);

    // Decompose: Jacobian blocks
    int offset = n_f;
    for (const int sz : jacobian_sizes) {
        vec_basic block(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);
        reduced_exprs_J.push_back(block);
        offset += sz;
    }
}

void SymbolicScalarsVectors::funcjachess2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,
     std::vector<vec_basic> &reduced_exprs_J, std::vector<vec_basic> &reduced_exprs_H,
     const std::vector<Expression> &f, const std::vector<std::vector<Expression>>& inputs_J,
     const std::vector<std::vector<Expression>>& inputs_H) {

    vec_basic exprs;
    int count_f = f.size();
    int count_J = 0;
    int count_H = 0;

    // Add original function expressions
    for (const auto &fi : f) {
        exprs.push_back(fi.get_basic());
    }

    // Compute and append all Jacobian entries for each input group
    std::vector<int> J_sizes;
    for (const auto& input_vec : inputs_J) {
        int m = f.size();
        int n = input_vec.size();
        J_sizes.push_back(m * n);
        count_J += m * n;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                exprs.push_back(f[i].diff(input_vec[j]).get_basic());
            }
        }
    }

    // Compute and append all Hessian entries for each input group
    std::vector<int> H_sizes;
    for (const auto& input_vec : inputs_H) {
        int m = f.size();
        int n = input_vec.size();
        H_sizes.push_back(m * n * n);
        count_H += m * n * n;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                Expression first = f[i].diff(input_vec[j]);
                for (int k = 0; k < n; ++k) {
                    exprs.push_back(first.diff(input_vec[k]).get_basic());
                }
            }
        }
    }

    // Apply Common Subexpression Elimination
    vec_basic reduced_exprs;
    cse(replacements, reduced_exprs, exprs);

    // Decompose reduced_exprs into f, J, and H
    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + count_f);

    int offset = count_f;
    for (int sz : J_sizes) {
        vec_basic Jblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);
        reduced_exprs_J.push_back(Jblock);
        offset += sz;
    }

    for (int sz : H_sizes) {
        vec_basic Hblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);
        reduced_exprs_H.push_back(Hblock);
        offset += sz;
    }
}

void SymbolicScalarsVectors::func2cppfiles(const std::vector<Expression> &f, const int functionid) {
    vec_pair replacements;
    vec_basic reduced_exprs;
    func2cse(replacements, reduced_exprs, f);

    // Determine variable usage
    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
    for (const auto &expr : f) {
        auto symbols = free_symbols(*expr.get_basic());
        used.insert(symbols.begin(), symbols.end());
    }

    auto depends_on = [&](const Expression &sym) {
        return used.count(sym.get_basic()) > 0;
    };


    // Generate function prototype header based on functionid
    std::ofstream hfile(funcnames[functionid] + std::string(".h"));
    hfile << "#pragma once\n\n";
    hfile << funcdecls[functionid] << ";\n";
    hfile.close();


    // Generate function prototype source based on functionid
    std::ofstream cppfile(funcnames[functionid] + std::string(".cpp"));
    cppfile << "#include \"" << funcnames[functionid] << ".h\"\n";
    cppfile << "#include <cmath>\n\n";
    cppfile << funcdecls[functionid] << "\n";
    cppfile << "{\n\n";

    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
    C99CodePrinter cpp;
    cppfile << "  for (int i = 0; i < N; ++i) {\n";
    
    // Emit symbolic variable loads
    for (const auto &[name, vec] : inputs) {
        for (size_t j = 0; j < vec.size(); ++j) {
            if (depends_on(vec[j]))
                cppfile << "    dstype " << name << j << " = " << name << "[" << j << "*N+i];\n";
        }
    }
    
    cppfile << "\n";
    
    // Emit intermediate CSE substitutions
    for (size_t n = 0; n < replacements.size(); ++n) {
        std::string var_name = cpp.apply(*replacements[n].first);
        std::string rhs = cpp.apply(*replacements[n].second);
        cppfile << "    dstype " << var_name << " = " << rhs << ";\n";
    }
    cppfile << "\n";
    
    for (size_t n = 0; n < f.size(); ++n) {
        cppfile << "    f[" << n << " * N + i] = " << cpp.apply(*reduced_exprs[n]) << ";\n";
    }
    
    cppfile << "  }\n";
    cppfile << "}\n\n";
    cppfile.close();
}

void SymbolicScalarsVectors::funcjac2cppfiles(const std::vector<Expression> &f, const int functionid) {
    vec_pair replacements;
    vec_basic reduced_exprs_f;
    std::vector<vec_basic> reduced_exprs_J;
    funcjac2cse(replacements, reduced_exprs_f, reduced_exprs_J, f, jacobianInputs[functionid]);

    // Determine variable usage
    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
    for (const auto &expr : f) {
        auto symbols = free_symbols(*expr.get_basic());
        used.insert(symbols.begin(), symbols.end());
    }

    auto depends_on = [&](const Expression &sym) {
        return used.count(sym.get_basic()) > 0;
    };

    // Generate function prototype header based on functionid
    std::ofstream hfile(funcnames[functionid] + std::string(".h"), std::ios::app);
    hfile << "void " << funcnames[functionid] << "jac" << "(dstype* f, ";
    int nJ = reduced_exprs_J.size();
    for (int k = 0; k < nJ; ++k)
        hfile << "dstype* J" << (k+1) << ", ";
    hfile << funcjacdecls[functionid] << ";\n";
    hfile.close();

    std::ofstream cppfile(funcnames[functionid] + std::string(".cpp"), std::ios::app);
    cppfile << "void " << funcnames[functionid] << "jac" << "(dstype* f, ";
    for (int k = 0; k < nJ; ++k)
        cppfile << "dstype* J" << (k+1) << ", ";
    cppfile << funcjacdecls[functionid] << "\n";
    cppfile << "{\n\n";

    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
    C99CodePrinter cpp;
    cppfile << "  for (int i = 0; i < N; ++i) {\n";

    for (const auto &[name, vec] : inputs) {
        for (size_t j = 0; j < vec.size(); ++j) {
            if (depends_on(vec[j]))
                cppfile << "    dstype " << name << j << " = " << name << "[" << j << "*N+i];\n";
        }
    }

    for (size_t n = 0; n < replacements.size(); ++n) {
        std::string var_name = cpp.apply(*replacements[n].first);
        std::string rhs = cpp.apply(*replacements[n].second);
        cppfile << "    dstype " << var_name << " = " << rhs << ";\n";
    }
    cppfile << "\n";
    
    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {
        cppfile << "    f[" << n << " * N + i] = " << cpp.apply(*reduced_exprs_f[n]) << ";\n";
    }

    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {
        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {
            cppfile << "    J" << (k+1) << "[" << j << " * N + i] = " << cpp.apply(*reduced_exprs_J[k][j]) << ";\n";
        }
    }
    cppfile << "  }\n";
    cppfile << "}\n\n";
    cppfile.close();
}

void SymbolicScalarsVectors::funcjachess2cppfiles(const std::vector<Expression> &f, const int functionid) {
    vec_pair replacements;
    vec_basic reduced_exprs_f;
    std::vector<vec_basic> reduced_exprs_J, reduced_exprs_H;
    funcjachess2cse(replacements, reduced_exprs_f, reduced_exprs_J, reduced_exprs_H, f, jacobianInputs[functionid], hessianInputs[functionid]);

    // Determine variable usage
    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
    for (const auto &expr : f) {
        auto symbols = free_symbols(*expr.get_basic());
        used.insert(symbols.begin(), symbols.end());
    }

    auto depends_on = [&](const Expression &sym) {
        return used.count(sym.get_basic()) > 0;
    };

    // Generate function prototype header based on functionid
    std::ofstream hfile(funcnames[functionid] + std::string(".h"), std::ios::app);
    hfile << "void " << funcnames[functionid] << "jachess" << "(dstype* f, ";
    int nJ = reduced_exprs_J.size();
    int nH = reduced_exprs_H.size();
    for (int k = 0; k < nJ; ++k) hfile << "dstype* J" << (k+1) << ", ";
    for (int k = 0; k < nH; ++k) hfile << "dstype* H" << (k+1) << ", ";
    hfile << funcjacdecls[functionid] << ";\n";
    hfile.close();

    std::ofstream cppfile(funcnames[functionid] + std::string(".cpp"), std::ios::app);
    cppfile << "void " << funcnames[functionid] << "jachess" << "(dstype* f, ";
    for (int k = 0; k < nJ; ++k) cppfile << "dstype* J" << (k+1) << ", ";
    for (int k = 0; k < nH; ++k) cppfile << "dstype* H" << (k+1) << ", ";
    cppfile << funcjacdecls[functionid] << "\n";
    cppfile << "{\n\n";

    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
    C99CodePrinter cpp;
    cppfile << "  for (int i = 0; i < N; ++i) {\n";

    for (const auto &[name, vec] : inputs) {
        for (size_t j = 0; j < vec.size(); ++j) {
            if (depends_on(vec[j]))
                cppfile << "    dstype " << name << j << " = " << name << "[" << j << "*N+i];\n";
        }
    }

    for (size_t n = 0; n < replacements.size(); ++n) {
        std::string var_name = cpp.apply(*replacements[n].first);
        std::string rhs = cpp.apply(*replacements[n].second);
        cppfile << "    dstype " << var_name << " = " << rhs << ";\n";
    }
    cppfile << "\n";
    
    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {
        cppfile << "    f[" << n << " * N + i] = " << cpp.apply(*reduced_exprs_f[n]) << ";\n";
    }

    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {
        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {
            cppfile << "    J" << (k+1) << "[" << j << " * N + i] = " << cpp.apply(*reduced_exprs_J[k][j]) << ";\n";
        }
    }

    for (size_t k = 0; k < reduced_exprs_H.size(); ++k) {
        for (size_t j = 0; j < reduced_exprs_H[k].size(); ++j) {
            cppfile << "    H" << (k+1) << "[" << j << " * N + i] = " << cpp.apply(*reduced_exprs_H[k][j]) << ";\n";
        }
    }
    cppfile << "  }\n";
    cppfile << "}\n\n";
    cppfile.close();
}

