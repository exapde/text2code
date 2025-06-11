#include <cstdlib>
#include <iostream>
#include <sstream>
#include <filesystem>
#include "CodeGenerator.cpp"

// g++ -O2 -std=c++17 Text2Code.cpp -o text2code

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./text2code <pdemodel.txt> [<symengine_path>]\n";
        return 1;
    }    
    
    try {
        if (std::filesystem::exists(argv[1])) {
            std::cout << "Generating the C++ code for this input text file ("<< argv[1] << ") ... \n\n";
        } else {
            std::cout << "Input file does not exist.\n";
            return 1;
        }  
        
        ParsedSpec spec = TextParser::parseFile(argv[1]);

        // Generate symbolic code
        CodeGenerator gen(spec);
        gen.generateCode2Cpp("Code2Cpp.cpp");
        gen.generateSymbolicFunctionsHpp("SymbolicFunctions.hpp");
        gen.generateSymbolicFunctionsCpp("SymbolicFunctions.cpp");        
        gen.generateSymbolicScalarsVectorsHpp("SymbolicScalarsVectors.hpp");
        gen.generateSymbolicScalarsVectorsCpp("SymbolicScalarsVectors.cpp");   
        
        std::string framework = spec.framework;
        if ((framework == "Cuda") || (framework == "CUDA") || (framework == "cuda")) {
          gen.generateCudaHipHpp("CudaHip.hpp");
        } 
        else if ((framework == "Hip") || (framework == "HIP") || (framework == "hip")) {
          gen.generateCudaHipHpp("CudaHip.hpp");
        } 

        std::cout << "The C++ code is generated in the following files:\n";
        std::cout << "    Code2Cpp.cpp\n";
        std::cout << "    SymbolicFunctions.hpp and SymbolicFunctions.cpp\n";
        std::cout << "    SymbolicScalarsVectors.hpp and SymbolicScalarsVectors.cpp\n";
        
        std::filesystem::path cwd = std::filesystem::current_path();
        std::cout << "These files are located in this directory: " << cwd << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error during parsing or generation: " << ex.what() << "\n";
        return 1;
    }

    // You should replace "/Users/ngoccuongnguyen/GitHub/text2code" with the correct symengine path on your system
    std::string symengine_path = (argc >= 3) ? argv[2] : "/Users/ngoccuongnguyen/GitHub/text2code";
    
    if (std::filesystem::exists(symengine_path) && std::filesystem::is_directory(symengine_path)) {
        std::cout << "\nCompiling the C++ generated code (Code2Cpp.cpp) ...\n";
    } else {
        std::cout <<"This symengine path ("<<symengine_path<<") is not valid. Please open Text2Code.cpp and replace it.\n";
        return 1;
    }
    
    // Construct compile command using symengine_path
    std::stringstream cmd;
    cmd << "g++ -std=c++17 -Wno-inconsistent-missing-override "
        << "Code2Cpp.cpp SymbolicFunctions.cpp SymbolicScalarsVectors.cpp "
        << "-I" << symengine_path << "/include "
        << "-I" << symengine_path << " "
        << symengine_path << "/lib/libsymengine.a "
        << "-o code2cpp";
    
    int status = std::system(cmd.str().c_str());
    if (status != 0) {
        std::cerr << "Compilation failed!\n";
        return 1;
    }

    std::cout << "\nRunning the C++ generated code (Code2Cpp) ...\n\n";
    status = std::system("./code2cpp");
    if (status != 0) {
        std::cerr << "Generator execution failed!\n";
        return 1;
    }
    else {
        std::filesystem::path cwd = std::filesystem::current_path();
        std::cout << "The below source codes are generated and located in this directory: " << cwd << "\n";      
        ParsedSpec spec = TextParser::parseFile(argv[1]);
        for (const auto& funcname : spec.outputs) {
          std::cout<<funcname<<".h, "<<funcname<<".cpp\n"; 
        }
    }

    return 0;
}

