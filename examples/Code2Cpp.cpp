#include "SymbolicScalarsVectors.hpp"

int main() 
{
  SymbolicScalarsVectors ssv;

  for (int i=0; i<ssv.outputfunctions.size(); i++) {
    if (ssv.outputfunctions[i] == true) {
      std::vector<Expression> f = ssv.evaluateSymbolicFunctions(i);
      ssv.func2cppfiles(f, i);
      if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(f, i);
      if (ssv.hessianInputs[i].size() > 0) ssv.funcjachess2cppfiles(f, i);
    }
  }
}
