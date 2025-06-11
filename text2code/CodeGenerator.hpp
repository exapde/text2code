#pragma once
#include "TextParser.hpp"

class CodeGenerator {
public:
    CodeGenerator(const ParsedSpec& spec);

    void generateCode2Cpp(const std::string& filename) const;
    void generateSymbolicFunctionsHpp(const std::string& filename) const;
    void generateSymbolicFunctionsCpp(const std::string& filename) const;
    void generateSymbolicScalarsVectorsHpp(const std::string& filename) const;
    void generateSymbolicScalarsVectorsCpp(const std::string& filename) const;
    void generateCudaHipHpp(const std::string& filename) const;
    
private:
    const ParsedSpec& spec;

    void generateFunctionHeader(std::ostream& os, const FunctionDef& func) const;
    void generateFunctionSource(std::ostream& os, const FunctionDef& func) const;
};




