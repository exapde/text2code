#include "CodeGenerator.hpp"

std::string prefixSymEngineFunctions(const std::string& expr) {
    static const std::regex func_regex(
    //R"(\b(sin|cos|tan|cot|asin|acos|atan|acot|sinh|cosh|tanh|coth|asinh|acosh|atanh|acoth|sech|csch|log|log10|exp|sqrt|cbrt|abs|ceil|floor|pow|max|min|sign|pi)\b)");
            R"(\b(pi)\b)");
    return std::regex_replace(expr, func_regex, "SymEngine::$1");
}

CodeGenerator::CodeGenerator(const ParsedSpec& spec_) : spec(spec_) {}

void CodeGenerator::generateCode2Cpp(const std::string& filename) const {
    std::ofstream os(filename);

    os << "#include \"SymbolicScalarsVectors.hpp\"\n\n";
    os << "int main() \n";
    os << "{\n";
    os << "  SymbolicScalarsVectors ssv;\n\n";
    os << "  for (int i=0; i<ssv.outputfunctions.size(); i++) {\n";
    os << "    if (ssv.outputfunctions[i] == true) {\n";
    os << "      std::vector<Expression> f = ssv.evaluateSymbolicFunctions(i);\n";
    os << "      ssv.func2cppfiles(f, i);\n";
    os << "      if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(f, i);\n";
    os << "      if (ssv.hessianInputs[i].size() > 0) ssv.funcjachess2cppfiles(f, i);\n";
    os << "    }\n";
    os << "  }\n";
    os << "}\n";
    
    os.close();  
}

void CodeGenerator::generateCudaHipHpp(const std::string& filename) const {  
  std::ofstream os(filename);
  
  os << "#pragma once\n\n";
  os << "#include <string>\n\n";
  os << "#if defined(__HIPCC__)\n";
  os << "  #include <hip/hip_runtime.h>\n";
  os << "  #define GPU_DEVICE __device__\n";
  os << "  #define GPU_GLOBAL __global__\n";
  os << "  #define gpuDeviceSynchronize hipDeviceSynchronize\n";
  os << "  #define GPU_BACKEND \"HIP\"\n";
  os << "#elif defined(__CUDACC__)\n";
  os << "  #include <cuda_runtime.h>\n";
  os << "  #define GPU_DEVICE __device__\n";
  os << "  #define GPU_GLOBAL __global__\n";
  os << "  #define gpuDeviceSynchronize cudaDeviceSynchronize\n";
  os << "  #define GPU_BACKEND \"CUDA\"\n";
  os << "#else\n";
  os << "  #error \"GPU backend not recognized. Use nvcc or hipcc.\"\n";
  os << "#endif\n\n";
  os << "#define GPU_LAMBDA [=] GPU_DEVICE\n\n";
  os << "template <typename Functor>\n";
  os << "GPU_GLOBAL void lambda_kernel(Functor f, size_t N) {\n";
  os << "    size_t i = blockIdx.x * blockDim.x + threadIdx.x;\n";
  os << "    if (i < N) f(i);\n";
  os << "}\n\n";
  os << "template <typename Functor>\n";
  os << "void parallel_for(const std::string& label, size_t N, Functor f, int blockSize = 256) {\n";
  os << "    int numBlocks = (N + blockSize - 1) / blockSize;\n";
  os << "    lambda_kernel<<<numBlocks, blockSize>>>(f, N);\n";
  os << "    gpuDeviceSynchronize();\n";
  os << "}\n";
  
  os.close();  
}

void CodeGenerator::generateSymbolicFunctionsHpp(const std::string& filename) const {
    std::ofstream os(filename);
    os << "#pragma once\n\n";
    
    os << "#include <vector>\n";
    os << "#include <string>\n";
    os << "#include <fstream>\n";
    os << "#include <iostream>\n";
    os << "#include <symengine/expression.h>\n";
    os << "#include <symengine/functions.h>\n";
    os << "#include <symengine/printers/codegen.h>\n";
    //os << "#include <symengine/cse.h>\n";
    os << "#include <symengine/symbol.h>\n";
    os << "#include <symengine/matrix.h>\n\n";

    os << "using SymEngine::Expression;\n";
    os << "using SymEngine::vec_pair;\n";
    os << "using SymEngine::vec_basic;\n";
    os << "using SymEngine::symbol;\n";
    os << "using SymEngine::RCP;\n";
    os << "using SymEngine::Basic;\n";
    os << "using SymEngine::map_basic_basic;\n";
    os << "using SymEngine::CodePrinter;\n";
    os << "using SymEngine::C99CodePrinter;\n";
    
    os << "\n";
    for (const auto& func : spec.functions) {
        generateFunctionHeader(os, func);
    }
    
    os.close();  
}

void CodeGenerator::generateSymbolicFunctionsCpp(const std::string& filename) const {
    std::ofstream os(filename);
    os << "#include \"SymbolicFunctions.hpp\"\n\n";
    
    for (const auto& func : spec.functions) {
        generateFunctionSource(os, func);
    }
    
    os.close();  
}

void CodeGenerator::generateFunctionHeader(std::ostream& os, const FunctionDef& func) const {
    os << "std::vector<Expression> " << func.name << "(";
    for (size_t i = 0; i < func.args.size(); ++i) {
        bool isscalar = false;       
        for (int j=0; j<spec.scalars.size(); j++) {               
          if (spec.scalars[j] == func.args[i]) {
            isscalar = true;
            break;
          }
        }
        bool isvector = false;       
        for (int j=0; j<spec.namevectors.size(); j++) {               
          if (spec.namevectors[j] == func.args[i]) {
            isvector = true;
            break;
          }
        }
        
        if ((isscalar==false) && (isvector==false)) {
          printf("Error: function argument %s is not listed in scalars or vectors\n", func.args[i].c_str());
          exit(-1);
        }
        
        if (isscalar==true) os << "const Expression& " << func.args[i];
        if (isvector==true) os << "const std::vector<Expression>& " << func.args[i];           
      
        if (i + 1 < func.args.size()) os << ", ";
    }
    os << ");\n";
}

void CodeGenerator::generateFunctionSource(std::ostream& os, const FunctionDef& func) const {
    std::string lhs = func.output;
    int lhs_size = func.outputsize;
    os << "std::vector<Expression> " << func.name << "(";
    for (size_t i = 0; i < func.args.size(); ++i) {
        bool isscalar = false;       
        for (int j=0; j<spec.scalars.size(); j++) {               
          if (spec.scalars[j] == func.args[i]) {
            isscalar = true;
            break;
          }
        }
        bool isvector = false;       
        for (int j=0; j<spec.namevectors.size(); j++) {               
          if (spec.namevectors[j] == func.args[i]) {
            isvector = true;
            break;
          }
        }
        
        if ((isscalar==false) && (isvector==false)) {
          printf("Error: function argument %s is not listed in scalars or vectors\n", func.args[i].c_str());
          exit(-1);
        }
        
        if (isscalar==true) os << "const Expression& " << func.args[i];
        if (isvector==true) os << "const std::vector<Expression>& " << func.args[i];           
        
        if (i + 1 < func.args.size()) os << ", ";
    }
    os << ") {\n";

    os << "    std::vector<Expression> " << lhs << ";\n";
    os << "    " << lhs << ".resize(" << lhs_size << ");\n\n";
    
    std::regex for_loop_start(R"(^\s*for\s+(\w+)\s+in\s+(\d+):(\d+)\s*$)");
    std::regex for_loop_end(R"(^\s*endfor\s*$)");
    bool in_loop = false;

    std::regex call_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\((.*)\)\s*$)");
    std::regex assign_pattern(R"(^\s*(\w+)\s*=\s*(.+)$)");    
    std::regex zeros_pattern(R"(^\s*zeros\s*\(\s*(\w+)\s*,\s*(\d+)\s*\)\s*$)");
    std::regex  ones_pattern(R"(^\s*ones\s*\(\s*(\w+)\s*,\s*(\d+)\s*\)\s*$)");
    //std::regex fill_pattern(R"(^\s*fill\((\w+),\s*(\d+)\s*,\s*(\d+)\)\s*;?\s*$)");
    std::regex fill_pattern(R"(^\s*fill\((\w+),\s*(\d+),\s*([-+]?[0-9]*\.?[0-9]+)\)\s*;?\s*$)");
    std::regex zeros_no_size_pattern(R"(^\s*zeros\((\w+)\)\s*;?\s*$)");
    std::regex ones_no_size_pattern(R"(^\s*ones\((\w+)\)\s*;?\s*$)");
    std::regex fill_no_size_pattern(R"(^\s*fill\((\w+),\s*([-+]?[0-9]*\.?[0-9]+)\)\s*;?\s*$)");

    std::regex mat_decl_pattern(R"(^\s*matrix\s+(\w+)\((\d+),(\d+)\)\s*$)");
    std::regex matset_pattern(R"((\w+)\[(\d+)\]\[(\d+)\]\s*=\s*(.+))");
    std::regex matget_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\[(\d+)\]\[(\d+)\]\s*$)");
    std::regex det_pattern(R"(^\s*(\w+)\[(\d+)\]\s*=\s*det\((\w+)\)\s*$)");
    std::regex trace_pattern(R"(^\s*(\w+)\[(\d+)\]\s*=\s*trace\((\w+)\)\s*$)");
    std::regex inv_pattern(R"(^\s*(\w+)\s*=\s*inv\((\w+)\)\s*$)");
    std::regex transpose_pattern(R"(^\s*(\w+)\s*=\s*transpose\((\w+)\)\s*$)");
    std::regex mat_binop_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\s*([\+\-])\s*(\w+)\s*$)");
    std::regex scalar_mat_mul_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\s*\*\s*(\w+)\s*$)");
    std::regex matmul_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\s*\*\s*(\w+)(\s*\*\s*(\w+))?\s*$)");

    for (const auto& line_raw : func.body) {
        std::string line = prefixSymEngineFunctions(line_raw);             
        
        line.erase(std::remove(line.begin(), line.end(), ';'), line.end());

        std::smatch match;

        if (std::regex_match(line, match, for_loop_start)) {
            std::string var = match[1];
            std::string start = match[2];
            std::string end = match[3];
            os << "    for (int " << var << " = " << start << "; " << var << " <= " << end << "; ++" << var << ") {\n";
            in_loop = true;
        } else if (std::regex_match(line, match, for_loop_end)) {
            os << "    }\n";
            in_loop = false;
        } else if (std::regex_match(line, match, mat_decl_pattern)) {
            std::string mat = match[1];
            std::string rows = match[2];
            std::string cols = match[3];
            os << "    SymEngine::DenseMatrix " << mat << "(" << rows << ", " << cols << ");\n";
        } else if (std::regex_match(line, match, matset_pattern)) {
            std::string mat = match[1];
            std::string i = match[2];
            std::string j = match[3];
            std::string rhs = match[4];
            os << "    " << mat << ".set(" << i << ", " << j << ", " << rhs << ");\n";
        } else if (std::regex_match(line, match, matget_pattern)) {
            std::string lhs = match[1];
            std::string mat = match[2];
            std::string i = match[3];
            std::string j = match[4];
            os << "    Expression " << lhs << " = " << mat << ".get(" << i << ", " << j << ");\n";
        } else if (std::regex_match(line, match, det_pattern)) {
            std::string vec = match[1];
            std::string index = match[2];
            std::string mat = match[3];
            os << "    " << vec << "[" << index << "] = " << mat << ".det();\n";
        } else if (std::regex_match(line, match, trace_pattern)) {
            std::string vec = match[1];
            std::string index = match[2];
            std::string mat = match[3];
            os << "    " << vec << "[" << index << "] = SymEngine::trace(" << mat << ");\n";
        } else if (std::regex_match(line, match, inv_pattern)) {
            std::string lhs_mat = match[1];
            std::string rhs_mat = match[2];
            if (func.matrices.count(lhs_mat) && func.matrices.count(rhs_mat)) {
                os << "    DenseMatrix " << lhs_mat << " = SymEngine::inv(" << rhs_mat << ");\n";
            }
        } else if (std::regex_match(line, match, transpose_pattern)) {
            std::string lhs_mat = match[1];
            std::string rhs_mat = match[2];
            if (func.matrices.count(lhs_mat) && func.matrices.count(rhs_mat)) {
                os << "    DenseMatrix " << lhs_mat << " = SymEngine::transpose(" << rhs_mat << ");\n";
            }
        } else if (std::regex_match(line, match, mat_binop_pattern)) {
            std::string lhs_mat = match[1];
            std::string mat1 = match[2];
            std::string op = match[3];
            std::string mat2 = match[4];
            if (func.matrices.count(lhs_mat) && func.matrices.count(mat1) && func.matrices.count(mat2)) {
                os << "    DenseMatrix " << lhs_mat << " = DenseMatrix(" << mat1 << ") " << op << " DenseMatrix(" << mat2 << ");\n";
            } else {
                os << "    Expression " << lhs_mat << " = DenseMatrix(" << mat1 << ") " << op << " DenseMatrix(" << mat2 << ");\n";
            }
        } else if (std::regex_match(line, match, scalar_mat_mul_pattern)) {
            std::string lhs = match[1];
            std::string scalar = match[2];
            std::string mat = match[3];
            if (func.matrices.count(mat)) {
                os << "    DenseMatrix " << lhs << " = DenseMatrix(" << mat << ") * " << scalar << ";\n";
            } else {
                os << "    Expression " << lhs << " = " << scalar << " * " << mat << ";\n";
            }
        } else if (std::regex_match(line, match, matmul_pattern)) {
            std::string lhs = match[1];
            std::string m1 = match[2];
            std::string m2 = match[3];
            std::string m3 = match[5];
            bool lhs_is_matrix = func.matrices.count(lhs);
            bool rhs_has_matrix = func.matrices.count(m1) || func.matrices.count(m2) || (!m3.empty() && func.matrices.count(m3));
            if (lhs_is_matrix || rhs_has_matrix) {
                os << "    DenseMatrix " << lhs << " = " << m1 << " * " << m2;
                if (!m3.empty()) os << " * " << m3;
                os << ";\n";
            } else {
                os << "    Expression " << lhs << " = " << m1 << " * " << m2;
                if (!m3.empty()) os << " * " << m3;
                os << ";\n";
            }
        } else if (std::regex_match(line, match, call_pattern)) {
            std::string var = match[1];
            std::string fname = match[2];
            std::string args = match[3];
            if (var != lhs)
              os << "    auto " << var << " = " << fname << "(" << args << ");\n";
            else
              os << "    " << var << " = " << fname << "(" << args << ");\n";
        } else if (std::regex_match(line, match, zeros_pattern)) {
            std::string var = match[1];
            int size = std::stoi(match[2]);
            if (var != lhs) {
              os << "    std::vector<Expression> " << var << "(" << size << ", Expression(0));\n";              
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(0);\n";              
              os << "    }\n";
            }
        } else if (std::regex_match(line, match, zeros_no_size_pattern)) {
            std::string var = match[1];
            if (var != lhs) {
              printf("Error: zeros(%s) is invalid", var.c_str());
              exit(-1);
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << lhs_size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(0);\n";              
              os << "    }\n";
            }            
         } else if (std::regex_match(line, match, ones_no_size_pattern)) {
            std::string var = match[1];
            if (var != lhs) {
              printf("Error: ones(%s) is invalid", var.c_str());
              exit(-1);
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << lhs_size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(1);\n";              
              os << "    }\n";
            }                
        } else if (std::regex_match(line, match, ones_pattern)) {
            std::string var = match[1];
            int size = std::stoi(match[2]);
            if (var != lhs) {
              os << "    std::vector<Expression> " << var << "(" << size << ", Expression(1));\n";              
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(1);\n";              
              os << "    }\n";
            }
        } else if (std::regex_match(line, match, fill_pattern)) {
            std::string var = match[1];
            int size = std::stoi(match[2]);
            double value = std::stod(match[3]);
            if (var != lhs) {
              os << "    std::vector<Expression> " << var << "(" << size << ", Expression("<< value << "));\n";              
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression("<< value <<");\n";              
              os << "    }\n";
            }            
        } else if (std::regex_match(line, match, fill_no_size_pattern)) {
            std::string var = match[1];
            double value = std::stod(match[2]);
            if (var != lhs) {
              printf("Error: fill(%s, %g) is invalid", var.c_str(), value);
              exit(-1);
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << lhs_size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression("<< value <<");\n";              
              os << "    }\n";
            }
        } else if (line.find(lhs + "[") == 0) {
            size_t eq = line.find('=');
            std::string lhs_expr = line.substr(0, eq);
            std::string rhs_expr = line.substr(eq + 1);
            os << "    " << lhs_expr << " = " << rhs_expr << ";\n";
        } else if (std::regex_match(line, match, assign_pattern)) {
            std::string var = match[1];
            std::string rhs = match[2];
            os << "    Expression " << var << " = " << rhs << ";\n";
        }
    }

    os << "    return " << lhs << ";\n";
    os << "}\n\n";
}

void CodeGenerator::generateSymbolicScalarsVectorsHpp(const std::string& filename) const {
    std::ofstream os(filename);
    if (!os) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    os << "#pragma once\n\n";
    os << "#include \"SymbolicFunctions.hpp\"\n\n";
    
    os << "class SymbolicScalarsVectors {\n\n";
    os << "public:\n\n";
      
    // Scalars
    os << "    // input symbolic scalars\n";
    for (const auto& s : spec.scalars) {
        os << "    Expression " << s << ";\n";
    }
    os << "\n";

    // Vectors
    os << "    // input symbolic vectors\n";
    for (const auto& [name, size] : spec.vectors) {
        os << "    std::vector<Expression> " << name << ";\n";
    }
    os << "\n";

    // Vector sizes
    os << "    // vector sizes\n";
    for (const auto& [name, size] : spec.vectors) {
        os << "    int sz" << name << ";\n";
    }
    os << "\n";
        
    os << "    std::vector<bool> outputfunctions;\n";
    os << "    std::vector<std::vector<std::string>> funcargs;\n";
    os << "    std::vector<std::vector<std::string>> funcargssizes;\n";
    os << "    std::vector<std::string> funcnames;\n";
    os << "    std::vector<std::string> funcdecls;\n";
    os << "    std::vector<std::string> funcjacdecls;\n\n";
    
    os << "    std::vector<std::vector<std::pair<std::string, std::vector<Expression>>>> inputvectors;\n";  
    os << "    std::vector<std::vector<std::pair<std::string, Expression>>> inputscalars;\n";  
    os << "    std::vector<std::vector<std::vector<Expression>>> jacobianInputs;\n";
    os << "    std::vector<std::vector<std::vector<Expression>>> hessianInputs;\n\n";

    os << "    SymbolicScalarsVectors();\n\n";
    
    os << "    std::vector<Expression> evaluateSymbolicFunctions(int call);\n\n";    
        
    os << "    void func2cse(vec_pair &replacements, vec_basic &reduced_exprs, const std::vector<Expression> &f);\n\n";
    
    os << "    void funcjac2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,\n";
    os << "                         std::vector<vec_basic> &reduced_exprs_J, const std::vector<Expression> &f,\n";
    os << "                         const std::vector<std::vector<Expression>>& inputs_J);\n\n";

    os << "    void funcjachess2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,\n";
    os << "                         std::vector<vec_basic> &reduced_exprs_J, std::vector<vec_basic> &reduced_exprs_H,\n";
    os << "                         const std::vector<Expression> &f, const std::vector<std::vector<Expression>>& inputs_J,\n";
    os << "                         const std::vector<std::vector<Expression>>& inputs_H);\n\n";
    
    os << "    void func2cppfiles(const std::vector<Expression> &f, const int functionid);\n";
    os << "    void funcjac2cppfiles(const std::vector<Expression> &f, const int functionid);\n";
    os << "    void funcjachess2cppfiles(const std::vector<Expression> &f, const int functionid);\n";
    
    os << "};\n";
    
    os.close();  
}

void emitSymbolicScalarsVectors(std::ostream& os, const ParsedSpec& spec) {
  
    os << "SymbolicScalarsVectors::SymbolicScalarsVectors() {\n\n";

    for (const auto& s : spec.scalars) {
        os <<"    "<<s<< " =   Expression(\"" << s << "\")" << ";\n";
    }
    os << "\n";
                
    for (const auto& vec : spec.vectors) {
        const std::string& name = vec.first;
        int size = vec.second;
        std::string szname = "sz" + name;

        os << "    " << szname << " = " << size << ";\n";
        os << "    " << name << ".resize(" << size << ");\n";
        
        std::string iv = "i";
        os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
        os << "         " << name << "[" << iv << "] = Expression(\"" << name <<"\"  + std::to_string(i));\n";
        os << "    }\n";
                      
        os << "\n";
    }

    os << "    outputfunctions.assign(" << spec.functions.size() << ", false);\n";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
        const std::string& funcname = spec.functions[i].name;
        if (std::find(spec.outputs.begin(), spec.outputs.end(), funcname) != spec.outputs.end()) {
            os << "    outputfunctions[" << i << "] = true;\n";
        }
    }
    os << "\n";

    // Emit funcnames
    os << "    funcnames = {";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
        os<<"\""<<spec.functions[i].name<<"\"";
        if (i < spec.functions.size() - 1) os << ", ";
    }
    os << "};\n\n";
    
    // Emit funcargs
    os << "    funcargs = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";
        for (size_t i = 0; i < func.args.size(); ++i) {
            os << "\"" << func.args[i] << "\"";
            if (i < func.args.size() - 1) os << ", ";
        }
        os << "}";
        if (fi < spec.functions.size() - 1) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    // Emit funcargs
    os << "    funcargssizes = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";
        for (size_t i = 0; i < func.args.size(); ++i) {
            os << "\"sz" << func.args[i] << "\"";
            if (i < func.args.size() - 1) os << ", ";
        }
        os << "}";
        if (fi < spec.functions.size() - 1) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    os << "    funcdecls = {\n";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
      const FunctionDef& func = spec.functions[i];
      const std::string& type = spec.datatype;
      const std::string& name = func.name;
            
      os <<"       \""<<"void " << name << "(" << type << "* f";
      // Function arguments
      for (const std::string& arg : func.args) {
          bool is_scalar = std::find(spec.scalars.begin(), spec.scalars.end(), arg) != spec.scalars.end();
          os << ", const " << type;
          if (!is_scalar) os << "*";
          os << " " << arg;
      }
      // Metadata arguments
      os << ", const int N, const int szf";
      for (const std::string& arg : func.args) {
          if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end()) {
              os << ", const int sz" << arg;
          }
      }
      os << ")"<<"\"";
      if (i < spec.functions.size() - 1) os << ", ";
      os << "\n";
    }
    os << "    };\n\n";
                
    os << "    funcjacdecls = {\n";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
      const FunctionDef& func = spec.functions[i];
      const std::string& type = spec.datatype;
      const std::string& name = func.name;
            
      os <<"       \"";
      // Function arguments
      for (int i=0; i < func.args.size(); i++) {
          std::string arg = func.args[i];
          bool is_scalar = std::find(spec.scalars.begin(), spec.scalars.end(), arg) != spec.scalars.end();
          if (i==0)
            os << "const " << type;
          else 
            os << ", const " << type;
          if (!is_scalar) os << "*";
          os << " " << arg;
      }
      // Metadata arguments
      os << ", const int N, const int szf";
      for (const std::string& arg : func.args) {
          if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end()) {
              os << ", const int sz" << arg;
          }
      }
      os << ")"<<"\"";
      if (i < spec.functions.size() - 1) os << ", ";
      os << "\n";
    }
    os << "    };\n\n";
    
    os << "    inputvectors = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";        

        // Collect vector arguments only (exclude scalars)
        std::vector<std::string> vector_args;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end()) {
                vector_args.push_back(arg);
            }
        }

        for (size_t i = 0; i < vector_args.size(); ++i) {
            const auto& name = vector_args[i];
            os << "{\"" << name << "\", " << name << "}";
            if (i + 1 < vector_args.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";

    os << "    inputscalars = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";        

        // Collect scalar arguments 
        std::vector<std::string> scalar_args;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) != spec.scalars.end()) {
                scalar_args.push_back(arg);
            }
        }

        for (size_t i = 0; i < scalar_args.size(); ++i) {
            const auto& name = scalar_args[i];
            os << "{\"" << name << "\", " << name << "}";
            if (i + 1 < scalar_args.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    os << "    jacobianInputs = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";

        // Collect vector arguments that are in the Jacobian list
        std::vector<std::string> jacvecs;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end() &&
                std::find(spec.jacobian.begin(), spec.jacobian.end(), arg) != spec.jacobian.end()) {
                jacvecs.push_back(arg);
            }
        }

        for (size_t i = 0; i < jacvecs.size(); ++i) {
            os << jacvecs[i];
            if (i + 1 < jacvecs.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    os << "    hessianInputs = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";

        // Collect vector arguments that are in the Jacobian list
        std::vector<std::string> jacvecs;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end() &&
                std::find(spec.hessian.begin(), spec.hessian.end(), arg) != spec.hessian.end()) {
                jacvecs.push_back(arg);
            }
        }

        for (size_t i = 0; i < jacvecs.size(); ++i) {
            os << jacvecs[i];
            if (i + 1 < jacvecs.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
            
    os << "}\n\n";
}

void emitevaluateSymbolicFunctions(std::ostream& os, const ParsedSpec& spec) {
    os << "std::vector<Expression> SymbolicScalarsVectors::evaluateSymbolicFunctions(int call)\n";
    os << "{\n";
    os << "  std::vector<Expression> f;\n\n";
    os << "  switch (call) {\n";

    for (size_t i = 0; i < spec.functions.size(); ++i) {
        const auto& sig = spec.functions[i];
        os << "    case " << i << ":\n";
        os << "      f = " << sig.name << "(";
        for (size_t j = 0; j < sig.args.size(); ++j) {
            os << sig.args[j];
            if (j < sig.args.size() - 1)
                os << ", ";
        }
        os << ");\n";
        os << "      break;\n";
    }

    os << "    default:\n";
    os << "      throw std::runtime_error(\"Invalid function call in evaluateSymbolicFunctions\");\n";
    os << "  }\n\n";
    os << "  return f;\n";
    os << "}\n\n";
}
            
void emitfunc2cse(std::ostream& os) {
    os<<"void SymbolicScalarsVectors::func2cse(vec_pair &replacements, vec_basic &reduced_exprs, const std::vector<Expression> &f) {\n\n"
        "   vec_basic exprs;\n"
        "   for (const auto &fi : f) {\n"
        "       exprs.push_back(fi.get_basic());\n"
        "   }\n"
        "   cse(replacements, reduced_exprs, exprs);\n"
        "}\n\n";
}

void emitfuncjac2cse(std::ostream& os) {
    os << "void SymbolicScalarsVectors::funcjac2cse(vec_pair &replacements,\n"
       << "                                         vec_basic &reduced_exprs_f,\n"
       << "                                         std::vector<vec_basic> &reduced_exprs_J,\n"
       << "                                         const std::vector<Expression> &f,\n"
       << "                                         const std::vector<std::vector<Expression>>& inputs_J) {\n"
       << "    vec_basic exprs;\n\n"
       << "    // Track original sizes\n"
       << "    int n_f = static_cast<int>(f.size());\n"
       << "    std::vector<int> jacobian_sizes;\n\n"
       << "    // Add original function expressions\n"
       << "    for (const auto &fi : f) {\n"
       << "        exprs.push_back(fi.get_basic());\n"
       << "    }\n\n"
       << "    // Append Jacobians and record sizes for each input group\n"
       << "    for (const auto& input_vec : inputs_J) {\n"
       << "        int m = f.size();\n"
       << "        int n = input_vec.size();\n"
       << "        jacobian_sizes.push_back(m * n);\n"
       << "        for (int i = 0; i < m; ++i) {\n"
       << "            for (int j = 0; j < n; ++j) {\n"
       << "                exprs.push_back(f[i].diff(input_vec[j]).get_basic());\n"
       << "            }\n"
       << "        }\n"
       << "    }\n\n"
       << "    // Apply CSE\n"
       << "    vec_basic reduced_exprs;\n"
       << "    cse(replacements, reduced_exprs, exprs);\n\n"
       << "    // Decompose: f\n"
       << "    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + n_f);\n\n"
       << "    // Decompose: Jacobian blocks\n"
       << "    int offset = n_f;\n"
       << "    for (const int sz : jacobian_sizes) {\n"
       << "        vec_basic block(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);\n"
       << "        reduced_exprs_J.push_back(block);\n"
       << "        offset += sz;\n"
       << "    }\n"
       << "}\n\n";
}

void emitfuncjachess2cse(std::ostream& os) {
    os << "void SymbolicScalarsVectors::funcjachess2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,\n"
       << "     std::vector<vec_basic> &reduced_exprs_J, std::vector<vec_basic> &reduced_exprs_H,\n"
       << "     const std::vector<Expression> &f, const std::vector<std::vector<Expression>>& inputs_J,\n"
       << "     const std::vector<std::vector<Expression>>& inputs_H) {\n\n"

       << "    vec_basic exprs;\n"
       << "    int count_f = f.size();\n"
       << "    int count_J = 0;\n"
       << "    int count_H = 0;\n\n"

       << "    // Add original function expressions\n"
       << "    for (const auto &fi : f) {\n"
       << "        exprs.push_back(fi.get_basic());\n"
       << "    }\n\n"

       << "    // Compute and append all Jacobian entries for each input group\n"
       << "    std::vector<int> J_sizes;\n"
       << "    for (const auto& input_vec : inputs_J) {\n"
       << "        int m = f.size();\n"
       << "        int n = input_vec.size();\n"
       << "        J_sizes.push_back(m * n);\n"
       << "        count_J += m * n;\n"
       << "        for (int i = 0; i < m; ++i) {\n"
       << "            for (int j = 0; j < n; ++j) {\n"
       << "                exprs.push_back(f[i].diff(input_vec[j]).get_basic());\n"
       << "            }\n"
       << "        }\n"
       << "    }\n\n"

       << "    // Compute and append all Hessian entries for each input group\n"
       << "    std::vector<int> H_sizes;\n"
       << "    for (const auto& input_vec : inputs_H) {\n"
       << "        int m = f.size();\n"
       << "        int n = input_vec.size();\n"
       << "        H_sizes.push_back(m * n * n);\n"
       << "        count_H += m * n * n;\n"
       << "        for (int i = 0; i < m; ++i) {\n"
       << "            for (int j = 0; j < n; ++j) {\n"
       << "                Expression first = f[i].diff(input_vec[j]);\n"
       << "                for (int k = 0; k < n; ++k) {\n"
       << "                    exprs.push_back(first.diff(input_vec[k]).get_basic());\n"
       << "                }\n"
       << "            }\n"
       << "        }\n"
       << "    }\n\n"

       << "    // Apply Common Subexpression Elimination\n"
       << "    vec_basic reduced_exprs;\n"
       << "    cse(replacements, reduced_exprs, exprs);\n\n"

       << "    // Decompose reduced_exprs into f, J, and H\n"
       << "    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + count_f);\n\n"

       << "    int offset = count_f;\n"
       << "    for (int sz : J_sizes) {\n"
       << "        vec_basic Jblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);\n"
       << "        reduced_exprs_J.push_back(Jblock);\n"
       << "        offset += sz;\n"
       << "    }\n\n"

       << "    for (int sz : H_sizes) {\n"
       << "        vec_basic Hblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);\n"
       << "        reduced_exprs_H.push_back(Hblock);\n"
       << "        offset += sz;\n"
       << "    }\n"
       << "}\n\n";
}

void forloopstart(std::ostream& os, std::string framework)
{       
    if ((framework == "kokkos") || (framework == "Kokkos") || (framework == "KOKKOS")) {    
      os << "    cppfile << \"  Kokkos::parallel_for(\\\"\"<< funcnames[functionid] <<\"\\\", N, KOKKOS_LAMBDA(const size_t i) {\\n\";\n";                  
    }        
    else if ((framework == "cuda") || (framework == "Cuda") || (framework == "CUDA")) {    
      os << "    cppfile << \"  parallel_for(\\\"\"<< funcnames[functionid] <<\"\\\", N, GPU_LAMBDA(const size_t i) {\\n\";\n";                  
    }        
    else if ((framework == "hip") || (framework == "Hip") || (framework == "HIP")) {    
      os << "    cppfile << \"  parallel_for(\\\"\"<< funcnames[functionid] <<\"\\\", N, GPU_LAMBDA(const size_t i) {\\n\";\n";                  
    }        
    else {
      os << "    cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n";   
    }
}

void forloopend(std::ostream& os, std::string framework)
{       
    if ((framework == "kokkos") || (framework == "Kokkos") || (framework == "KOKKOS")) {    
      os << "    cppfile << \"  });\\n\";\n";
    }
    else if ((framework == "cuda") || (framework == "Cuda") || (framework == "CUDA")) {       
      os << "    cppfile << \"  });\\n\";\n";
    }
    else if ((framework == "hip") || (framework == "Hip") || (framework == "HIP")) {        
      os << "    cppfile << \"  });\\n\";\n";
    }
    else {
      os << "    cppfile << \"  }\\n\";\n";   
    }
}

void emitCudaHipHpp(std::ostream& os, std::string framework)
{
    if ((framework == "cuda") || (framework == "Cuda") || (framework == "CUDA")) {             
      os << "    hfile << \"#include \\\"CudaHip.hpp\\\" \\n\\n\";\n";
    }
    else if ((framework == "hip") || (framework == "Hip") || (framework == "HIP")) {        
      os << "    hfile << \"#include \\\"CudaHip.hpp\\\" \\n\\n\";\n";
    }  
}

void emitfunc2cppfiles(std::ostream& os, const ParsedSpec& spec) {
  
    os << "void SymbolicScalarsVectors::func2cppfiles(const std::vector<Expression> &f, const int functionid) {\n";
    
    os << "    vec_pair replacements;\n";
    os << "    vec_basic reduced_exprs;\n";
    os << "    func2cse(replacements, reduced_exprs, f);\n\n";

    os << "    // Determine variable usage\n";
    os << "    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;\n";
    os << "    for (const auto &expr : f) {\n";
    os << "        auto symbols = free_symbols(*expr.get_basic());\n";
    os << "        used.insert(symbols.begin(), symbols.end());\n";
    os << "    }\n\n";
    
    os << "    auto depends_on = [&](const Expression &sym) {\n";
    os << "        return used.count(sym.get_basic()) > 0;\n";
    os << "    };\n\n";
    
    os << "\n";
    os << "    // Generate function prototype header based on functionid\n";
    os << "    std::ofstream hfile(funcnames[functionid] + std::string(\".h\"));\n";
    os << "    hfile << \"#pragma once\\n\\n\";\n";
    emitCudaHipHpp(os, spec.framework);
    os << "    hfile << funcdecls[functionid] << \";\\n\";\n";
    os << "    hfile.close();\n";
    os << "\n\n";
    
    os << "    // Generate function prototype source based on functionid\n";
    os << "    std::ofstream cppfile(funcnames[functionid] + std::string(\".cpp\"));\n";
    os << "    cppfile << \"#include \\\"\" << funcnames[functionid] << \".h\\\"\\n\";\n";
    os << "    cppfile << \"#include <cmath>\\n\\n\";\n";
    os << "    cppfile << funcdecls[functionid] << \"\\n\";\n";    
    os << "    cppfile << \"{\\n\\n\";";
    os << "\n\n";        

    os << "    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];\n";
    os << "    C99CodePrinter cpp;\n";
    
    //os << "    cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n";
    forloopstart(os, spec.framework);    
    os << "    \n";
    os << "    // Emit symbolic variable loads\n";
    os << "    for (const auto &[name, vec] : inputs) {\n";
    os << "        for (size_t j = 0; j < vec.size(); ++j) {\n";
    os << "            if (depends_on(vec[j]))\n";
    os << "                cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*N+i];\\n\";\n";
    os << "        }\n";
    os << "    }\n";
    os << "    \n";
    os << "    cppfile << \"\\n\";\n";
    os << "    \n";
    
//     for (size_t i = 0; i < sf.substitutions.size(); ++i) {
//         std::string var_name = cpp.apply(*sf.substitutions[i].first);
//         std::string rhs = cpp.apply(*sf.substitutions[i].second);
//         cppfile << "    double " << var_name << " = " << rhs << ";\n";
//     }
// 
                
    os << "    // Emit intermediate CSE substitutions\n";
    os << "    for (size_t n = 0; n < replacements.size(); ++n) {\n";
    os << "        std::string var_name = cpp.apply(*replacements[n].first);\n";
    os << "        std::string rhs = cpp.apply(*replacements[n].second);\n";
    os << "        cppfile << \"    " << spec.datatype << " \" << var_name << \" = \" << rhs << \";\\n\";\n";
    //os << "        cppfile << \"    " << spec.datatype << " tmp\" << n << \" = \" << rhs << \";\\n\";\n";
    os << "    }\n";
    os << "    cppfile << \"\\n\";\n";
    os << "    \n";
    
    os << "    for (size_t n = 0; n < f.size(); ++n) {\n";
    //os << "        auto replaced = reduced_exprs[n]->subs(rename_map);\n";
    os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*reduced_exprs[n]) << \";\\n\";\n";
    os << "    }\n";
    os << "    \n";
    //os << "    cppfile << \"  }\\n\";\n";  // forloop
    forloopend(os, spec.framework);
    
    os << "    cppfile << \"}\\n\\n\";\n"; // function 
    os << "    cppfile.close();\n";
    os << "}\n\n";
    
//     for (int i = 0; i < k; ++i) {
//         cppfile << "    fb[" << i << "] = " << cpp.apply(sf.f[i].get_basic()) << ";\n";
//         for (int j = 0; j < m; ++j)
//             cppfile << "    dfdx_b[" << i*m + j << "] = " << cpp.apply(sf.dfdx[i][j].get_basic()) << ";\n";
//         for (int j = 0; j < n; ++j)
//             cppfile << "    dfdy_b[" << i*n + j << "] = " << cpp.apply(sf.dfdy[i][j].get_basic()) << ";\n";
//     }
    
//     os << "    // Emit final expressions\n";
//     os << "    map_basic_basic rename_map;\n";
//     os << "    for (size_t n = 0; n < replacements.size(); ++n)\n";
//     os << "        rename_map[replacements[n].first] = symbol(\"tmp\" + std::to_string(n));\n";
//     os << "    \n";
//     os << "    for (size_t n = 0; n < f.size(); ++n) {\n";
//     os << "        auto replaced = reduced_exprs[n]->subs(rename_map);\n";
//     os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "    }\n";
//     os << "    \n";
//     os << "    cppfile << \"  }\\n\";\n";
//     os << "    cppfile << \"}\\n\\n\";\n";    
//     os << "}\n\n";
}

void emitfuncjac2cppfiles(std::ostream& os, const ParsedSpec& spec) {
    os << "void SymbolicScalarsVectors::funcjac2cppfiles(const std::vector<Expression> &f, const int functionid) {\n";
    os << "    vec_pair replacements;\n";
    os << "    vec_basic reduced_exprs_f;\n";
    os << "    std::vector<vec_basic> reduced_exprs_J;\n";
    os << "    funcjac2cse(replacements, reduced_exprs_f, reduced_exprs_J, f, jacobianInputs[functionid]);\n\n";
    
    os << "    // Determine variable usage\n";
    os << "    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;\n";
    os << "    for (const auto &expr : f) {\n";
    os << "        auto symbols = free_symbols(*expr.get_basic());\n";
    os << "        used.insert(symbols.begin(), symbols.end());\n";
    os << "    }\n\n";
    
    os << "    auto depends_on = [&](const Expression &sym) {\n";
    os << "        return used.count(sym.get_basic()) > 0;\n";
    os << "    };\n\n";

    os << "    // Generate function prototype header based on functionid\n";
    os << "    std::ofstream hfile(funcnames[functionid] + std::string(\".h\"), std::ios::app);\n";
    os << "    hfile << \"void \" << funcnames[functionid] << \"jac\" << \"(" << spec.datatype << "* f, \";\n";
    os << "    int nJ = reduced_exprs_J.size();\n";
    os << "    for (int k = 0; k < nJ; ++k)\n";
    os << "        hfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
    os << "    hfile << funcjacdecls[functionid] << \";\\n\";\n";    
    os << "    hfile.close();\n";
    os << "\n";    
        
    os << "    std::ofstream cppfile(funcnames[functionid] + std::string(\".cpp\"), std::ios::app);\n";
    
    os << "    cppfile << \"void \" << funcnames[functionid] << \"jac\" << \"(" << spec.datatype << "* f, \";\n";
    os << "    for (int k = 0; k < nJ; ++k)\n";
    os << "        cppfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
    os << "    cppfile << funcjacdecls[functionid] << \"\\n\";\n";    
    os << "    cppfile << \"{\\n\\n\";";
    os << "\n\n";
    
//     os << "    const std::vector<std::string>& args = funcargs[functionid];\n";
//     os << "    const std::vector<std::string>& sizes = funcargssizes[functionid];\n";
//     os << "    for (size_t i = 0; i < args.size(); ++i)\n";
//     os << "        cppfile << \"const double* \" << args[i] << \", \";\n";
//     os << "    cppfile << \"const int N, const int szf, \";\n";
//     os << "    for (size_t i = 0; i < sizes.size(); ++i) {\n";
//     os << "        cppfile << \"const int \" << sizes[i];\n";
//     os << "        if (i + 1 < sizes.size()) cppfile << \", \";\n";
//     os << "    }\n";
//     os << "    cppfile << \")\\n\";\n";
//     os << "    cppfile << \"{\\n\\n\";\n\n";
    
    os << "    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];\n";
    os << "    C99CodePrinter cpp;\n";    
    
    //os << "    cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n\n";
    forloopstart(os, spec.framework);    
    os << "    for (const auto &[name, vec] : inputs) {\n";
    os << "        for (size_t j = 0; j < vec.size(); ++j) {\n";
    os << "            if (depends_on(vec[j]))\n";
    os << "                cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*N+i];\\n\";\n";
    os << "        }\n";
    os << "    }\n\n";

//     os << "    for (size_t n = 0; n < replacements.size(); ++n) {\n";
//     os << "        std::string rhs = cpp.apply(*replacements[n].second);\n";
//     os << "        cppfile << \"    " << spec.datatype << " tmp\" << n << \" = \" << rhs << \";\\n\";\n";
//     os << "    }\n\n";
// 
//     os << "    map_basic_basic rename_map;\n";
//     os << "    for (size_t n = 0; n < replacements.size(); ++n)\n";
//     os << "        rename_map[replacements[n].first] = symbol(\"tmp\" + std::to_string(n));\n\n";
// 
//     os << "    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {\n";
//     os << "        auto replaced = reduced_exprs_f[n]->subs(rename_map);\n";
//     os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "    }\n\n";
// 
//     os << "    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {\n";
//     os << "        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {\n";
//     os << "            auto replaced = reduced_exprs_J[k][j]->subs(rename_map);\n";
//     os << "            cppfile << \"    J\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "        }\n";
//     os << "    }\n";

    os << "    for (size_t n = 0; n < replacements.size(); ++n) {\n";
    os << "        std::string var_name = cpp.apply(*replacements[n].first);\n";
    os << "        std::string rhs = cpp.apply(*replacements[n].second);\n";
    os << "        cppfile << \"    " << spec.datatype << " \" << var_name << \" = \" << rhs << \";\\n\";\n";
    os << "    }\n";
    os << "    cppfile << \"\\n\";\n";
    os << "    \n";

    // No rename_map or substitution
    os << "    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {\n";
    os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*reduced_exprs_f[n]) << \";\\n\";\n";
    os << "    }\n\n";

    os << "    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {\n";
    os << "        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {\n";
    os << "            cppfile << \"    J\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*reduced_exprs_J[k][j]) << \";\\n\";\n";
    os << "        }\n";
    os << "    }\n";    
    //os << "    cppfile << \"  }\\n\";\n";
    forloopend(os, spec.framework);
    
    os << "    cppfile << \"}\\n\\n\";\n";
    os << "    cppfile.close();\n";
    os << "}\n\n";
}

void emitfuncjachess2cppfiles(std::ostream& os, const ParsedSpec& spec) {
    os << "void SymbolicScalarsVectors::funcjachess2cppfiles(const std::vector<Expression> &f, const int functionid) {\n";
    os << "    vec_pair replacements;\n";
    os << "    vec_basic reduced_exprs_f;\n";
    os << "    std::vector<vec_basic> reduced_exprs_J, reduced_exprs_H;\n";
    os << "    funcjachess2cse(replacements, reduced_exprs_f, reduced_exprs_J, reduced_exprs_H, f, jacobianInputs[functionid], hessianInputs[functionid]);\n\n";

    os << "    // Determine variable usage\n";
    os << "    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;\n";
    os << "    for (const auto &expr : f) {\n";
    os << "        auto symbols = free_symbols(*expr.get_basic());\n";
    os << "        used.insert(symbols.begin(), symbols.end());\n";
    os << "    }\n\n";
    
    os << "    auto depends_on = [&](const Expression &sym) {\n";
    os << "        return used.count(sym.get_basic()) > 0;\n";
    os << "    };\n\n";

    os << "    // Generate function prototype header based on functionid\n";
    os << "    std::ofstream hfile(funcnames[functionid] + std::string(\".h\"), std::ios::app);\n";    
    os << "    hfile << \"void \" << funcnames[functionid] << \"jachess\" << \"(" << spec.datatype << "* f, \";\n";
    os << "    int nJ = reduced_exprs_J.size();\n";
    os << "    int nH = reduced_exprs_H.size();\n";
    os << "    for (int k = 0; k < nJ; ++k) hfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
    os << "    for (int k = 0; k < nH; ++k) hfile << \""<< spec.datatype <<"* H\" << (k+1) << \", \";\n";    
    os << "    hfile << funcjacdecls[functionid] << \";\\n\";\n";    
    os << "    hfile.close();\n";
    os << "\n";    
            
    os << "    std::ofstream cppfile(funcnames[functionid] + std::string(\".cpp\"), std::ios::app);\n";
    os << "    cppfile << \"void \" << funcnames[functionid] << \"jachess\" << \"(" << spec.datatype << "* f, \";\n";    
    os << "    for (int k = 0; k < nJ; ++k) cppfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
    os << "    for (int k = 0; k < nH; ++k) cppfile << \""<< spec.datatype <<"* H\" << (k+1) << \", \";\n";    
    os << "    cppfile << funcjacdecls[functionid] << \"\\n\";\n";    
    os << "    cppfile << \"{\\n\\n\";";
    os << "\n\n";

//     os << "    const std::vector<std::string>& args = funcargs[functionid];\n";
//     os << "    const std::vector<std::string>& sizes = funcargssizes[functionid];\n";
//     os << "    for (size_t i = 0; i < args.size(); ++i)\n";
//     os << "        cppfile << \"const double* \" << args[i] << \", \";\n";
//     os << "    cppfile << \"const int N, const int szf, \";\n";
//     os << "    for (size_t i = 0; i < sizes.size(); ++i) {\n";
//     os << "        cppfile << \"const int \" << sizes[i];\n";
//     os << "        if (i + 1 < sizes.size()) cppfile << \", \";\n";
//     os << "    }\n";
//     os << "    cppfile << \")\\n\";\n";
//     os << "    cppfile << \"{\\n\\n\";\n\n";

    os << "    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];\n";
    os << "    C99CodePrinter cpp;\n";    
    
    //os << "    cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n\n";
    forloopstart(os, spec.framework);    
    os << "    for (const auto &[name, vec] : inputs) {\n";
    os << "        for (size_t j = 0; j < vec.size(); ++j) {\n";
    os << "            if (depends_on(vec[j]))\n";
    os << "                cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*N+i];\\n\";\n";
    os << "        }\n";
    os << "    }\n\n";

//     os << "    for (size_t n = 0; n < replacements.size(); ++n) {\n";
//     os << "        std::string rhs = cpp.apply(*replacements[n].second);\n";
//     os << "        cppfile << \"    " << spec.datatype << " tmp\" << n << \" = \" << rhs << \";\\n\";\n";
//     os << "    }\n\n";
// 
//     os << "    map_basic_basic rename_map;\n";
//     os << "    for (size_t n = 0; n < replacements.size(); ++n)\n";
//     os << "        rename_map[replacements[n].first] = symbol(\"tmp\" + std::to_string(n));\n\n";
// 
//     os << "    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {\n";
//     os << "        auto replaced = reduced_exprs_f[n]->subs(rename_map);\n";
//     os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "    }\n\n";
// 
//     os << "    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {\n";
//     os << "        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {\n";
//     os << "            auto replaced = reduced_exprs_J[k][j]->subs(rename_map);\n";
//     os << "            cppfile << \"    J\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "        }\n";
//     os << "    }\n\n";
// 
//     os << "    for (size_t k = 0; k < reduced_exprs_H.size(); ++k) {\n";
//     os << "        for (size_t j = 0; j < reduced_exprs_H[k].size(); ++j) {\n";
//     os << "            auto replaced = reduced_exprs_H[k][j]->subs(rename_map);\n";
//     os << "            cppfile << \"    H\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "        }\n";
//     os << "    }\n";

    os << "    for (size_t n = 0; n < replacements.size(); ++n) {\n";
    os << "        std::string var_name = cpp.apply(*replacements[n].first);\n";
    os << "        std::string rhs = cpp.apply(*replacements[n].second);\n";
    os << "        cppfile << \"    " << spec.datatype << " \" << var_name << \" = \" << rhs << \";\\n\";\n";
    os << "    }\n";
    os << "    cppfile << \"\\n\";\n";
    os << "    \n";

    // Skip substitution, directly print the reduced expressions
    os << "    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {\n";
    os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*reduced_exprs_f[n]) << \";\\n\";\n";
    os << "    }\n\n";

    os << "    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {\n";
    os << "        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {\n";
    os << "            cppfile << \"    J\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*reduced_exprs_J[k][j]) << \";\\n\";\n";
    os << "        }\n";
    os << "    }\n\n";

    os << "    for (size_t k = 0; k < reduced_exprs_H.size(); ++k) {\n";
    os << "        for (size_t j = 0; j < reduced_exprs_H[k].size(); ++j) {\n";
    os << "            cppfile << \"    H\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*reduced_exprs_H[k][j]) << \";\\n\";\n";
    os << "        }\n";
    os << "    }\n";    
    //os << "    cppfile << \"  }\\n\";\n";
    forloopend(os, spec.framework);
    
    os << "    cppfile << \"}\\n\\n\";\n";
    os << "    cppfile.close();\n";
    os << "}\n\n";
}

void CodeGenerator::generateSymbolicScalarsVectorsCpp(const std::string& filename) const {
    std::ofstream os(filename);

    os << "#include \"SymbolicScalarsVectors.hpp\"\n\n";
        
    emitSymbolicScalarsVectors(os, spec);
    
    emitevaluateSymbolicFunctions(os, spec);
                           
    emitfunc2cse(os);
        
    emitfuncjac2cse(os);
    
    emitfuncjachess2cse(os);    
       
    emitfunc2cppfiles(os, spec);    
    
    emitfuncjac2cppfiles(os, spec);       
    
    emitfuncjachess2cppfiles(os, spec);
    
    os.close();  
}

