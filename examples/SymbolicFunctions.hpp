#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/printers/codegen.h>
#include <symengine/symbol.h>
#include <symengine/matrix.h>

using SymEngine::Expression;
using SymEngine::vec_pair;
using SymEngine::vec_basic;
using SymEngine::symbol;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::map_basic_basic;
using SymEngine::CodePrinter;
using SymEngine::C99CodePrinter;

std::vector<Expression> conductivity(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta);
std::vector<Expression> flux(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta);
std::vector<Expression> forloop(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta);
std::vector<Expression> mass(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta);
std::vector<Expression> source(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta);
std::vector<Expression> fbouhdg(const std::vector<Expression>& uq, const std::vector<Expression>& w, const std::vector<Expression>& v, const std::vector<Expression>& x, const Expression& time, const std::vector<Expression>& mu, const std::vector<Expression>& eta, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau);
std::vector<Expression> initu(const std::vector<Expression>& x);
