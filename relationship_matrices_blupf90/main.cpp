#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>
#include <set>
#include <sstream>
#include "matrix.h"

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <genotype_file> [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --ids <file>           ID subset file (one ID per line)" << std::endl;
    std::cout << "  --missing <value>      Missing value (default: 5)" << std::endl;
    std::cout << "  --matrices <list>      Comma-separated list of matrices to compute:" << std::endl;
    std::cout << "                         additive, dominance, additive_additive," << std::endl;
    std::cout << "                         dominance_dominance, additive_dominance, all" << std::endl;
    std::cout << "                         (default: all)" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << program_name << " geno.txt --matrices additive" << std::endl;
    std::cout << "  " << program_name << " geno.txt --matrices additive,dominance" << std::endl;
    std::cout << "  " << program_name << " geno.txt --ids subset.txt --missing 9 --matrices all" << std::endl;
    std::cout << std::endl;
    std::cout << "Note: When using --ids, the final matrix order will match the order in the IDs file." << std::endl;
}

std::set<std::string> parse_matrices(const std::string& matrices_str) {
    std::set<std::string> matrices;
    std::stringstream ss(matrices_str);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        // Remove espaços em branco
        item.erase(std::remove_if(item.begin(), item.end(), ::isspace), item.end());
        
        if (item == "all") {
            matrices.insert("additive");
            matrices.insert("dominance");
            matrices.insert("additive_additive");
            matrices.insert("dominance_dominance");
            matrices.insert("additive_dominance");
            break;
        } else if (item == "additive" || item == "dominance" || 
                   item == "additive_additive" || item == "dominance_dominance" || 
                   item == "additive_dominance") {
            matrices.insert(item);
        } else {
            std::cerr << "Warning: Unknown matrix type '" << item << "' ignored." << std::endl;
        }
    }
    
    return matrices;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string genofile = argv[1];
    std::string idsfile = "";
    int missingValue = 5;
    std::set<std::string> matrices_to_compute;
    
    // Parse argumentos
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--ids" && i + 1 < argc) {
            idsfile = argv[++i];
        } else if (arg == "--missing" && i + 1 < argc) {
            missingValue = std::stoi(argv[++i]);
        } else if (arg == "--matrices" && i + 1 < argc) {
            matrices_to_compute = parse_matrices(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    // Se nenhuma matriz foi especificada, calcular todas
    if (matrices_to_compute.empty()) {
        matrices_to_compute = {"additive", "dominance", "additive_additive", 
                              "dominance_dominance", "additive_dominance"};
    }
    
    try {
        // Cria o objeto Functions
        Functions func(genofile, idsfile, missingValue);
        
        std::cout << "File loaded successfully: " << genofile << std::endl;
        if (!idsfile.empty()) {
            std::cout << "ID subset file used: " << idsfile << std::endl;
            std::cout << "Matrix order will follow the ID file order." << std::endl;
        }
        std::cout << "Missing value set to: " << missingValue << std::endl;
        std::cout << "Dimensions: " << func.get_num_samples() << " samples x "
                  << func.get_num_markers() << " markers" << std::endl;
        
        // Mostrar os primeiros IDs na ordem final
        std::cout << "\nFinal ID order (first " << std::min(func.get_num_samples(), size_t(10)) << " IDs):" << std::endl;
        const auto &ids = func.get_ids();
        for (size_t i = 0; i < std::min(ids.size(), size_t(10)); i++) {
            std::cout << "  [" << i+1 << "] " << ids[i] << std::endl;
        }
        if (ids.size() > 10) {
            std::cout << "  ... (and " << (ids.size() - 10) << " more)" << std::endl;
        }
        
        std::cout << "\nComputing matrices:" << std::endl;
        for (const auto& matrix : matrices_to_compute) {
            std::cout << "  - " << matrix << std::endl;
        }
        std::cout << std::endl;
        
        // Calcular apenas as matrizes solicitadas
        if (matrices_to_compute.count("additive")) {
            std::cout << "Computing additive inverse matrix..." << std::endl;
            auto add_inv = func.get_inverse_additive_matrix(true, 0.01, 1e-8);
            func.save_vectorized(add_inv, "gi.txt");
            std::cout << "  Saved to: gi.txt (ordered as specified)" << std::endl;
        }
        
        if (matrices_to_compute.count("dominance")) {
            std::cout << "Computing dominance inverse matrix..." << std::endl;
            auto dom_inv = func.get_inverse_dominance_matrix(true, 0.01, 1e-8);
            func.save_vectorized(dom_inv, "di.txt");
            std::cout << "  Saved to: di.txt (ordered as specified)" << std::endl;
        }
        
        if (matrices_to_compute.count("additive_additive")) {
            std::cout << "Computing additive-additive inverse matrix..." << std::endl;
            auto add_add_inv = func.get_inverse_additive_additive_matrix(true, 0.01, 1e-8);
            func.save_vectorized(add_add_inv, "ggi.txt");
            std::cout << "  Saved to: ggi.txt (ordered as specified)" << std::endl;
        }
        
        if (matrices_to_compute.count("dominance_dominance")) {
            std::cout << "Computing dominance-dominance inverse matrix..." << std::endl;
            auto dom_dom_inv = func.get_inverse_dominance_dominance_matrix(true, 0.01, 1e-8);
            func.save_vectorized(dom_dom_inv, "ddi.txt");
            std::cout << "  Saved to: ddi.txt (ordered as specified)" << std::endl;
        }
        
        if (matrices_to_compute.count("additive_dominance")) {
            std::cout << "Computing additive-dominance inverse matrix..." << std::endl;
            auto add_dom_inv = func.get_inverse_additive_dominance_matrix(true, 0.01, 1e-8);
            func.save_vectorized(add_dom_inv, "gdi.txt");
            std::cout << "  Saved to: gdi.txt (ordered as specified)" << std::endl;
        }
        
        std::cout << "\nAll requested matrices computed successfully!" << std::endl;
        std::cout << "Matrix elements are indexed according to the final ID order shown above." << std::endl;
        
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
