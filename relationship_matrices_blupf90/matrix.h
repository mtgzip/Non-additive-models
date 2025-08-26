#ifndef MA_FUNCTIONS_H
#define MA_FUNCTIONS_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <cmath>

class Functions
{
private:
    // **LOAD FILES**
    std::string genofile; // File with genotype
    std::string idsfile;  // File with IDs for subset
    int MissingValue;     // Missing value in data

    // ** MAIN VECTORS **
    std::vector<std::vector<double>> snps; // Genotype data

    // ** DIMENSIONS **
    size_t num_samples;
    size_t num_markers;

    // ** RELATIONSHIP MATRIXS** //
    std::vector<std::vector<double>> additive_matrix;
    std::vector<std::vector<double>> dominance_matrix;
    std::vector<std::vector<double>> identity_matrix;
    std::vector<std::vector<double>> additive_matrix_inv;
    std::vector<std::vector<double>> dominance_matrix_inv;
    std::vector<std::vector<double>> identity_matrix_inv;
    
    // ** IDs **
    std::vector<std::string> ids_list;

    // ** SUBSET FUNCTIONALITY **
    std::vector<std::string> ordered_subset_ids; // Mantém a ordem dos IDs do arquivo
    std::unordered_set<std::string> subset_ids;  // Set of IDs to keep (para busca rápida)
    bool use_subset; // Flag to indicate if subset is being used

    // ** MATRIX OPERATIONS **
    std::vector<std::vector<double>> cholesky_decomposition(const std::vector<std::vector<double>> &matrix);
    std::vector<std::vector<double>> forward_substitution(const std::vector<std::vector<double>> &L, const std::vector<std::vector<double>> &B);
    std::vector<std::vector<double>> backward_substitution(const std::vector<std::vector<double>> &U, const std::vector<std::vector<double>> &B);
    std::vector<std::vector<double>> matrix_multiply(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);
    std::vector<std::vector<double>> matrix_transpose(const std::vector<std::vector<double>> &matrix);

    // ** HELPER FUNCTIONS **
    void load_subset_ids(const std::string &ids_file);

public:
    // Method declarations
    std::vector<std::vector<double>> loadgeno(const std::string &filename, int missingValue);
    std::pair<size_t, size_t> get_dimensions(const std::string &filename);

    // New function to process large genotype files with optional subset
    void process_large_genotypes(const std::string &file_path, const std::string &ids_file, int missingValue);

    // Constructor modificado para aceitar arquivo de IDs opcional
    Functions(const std::string &_genofile, const std::string &_idsfile, int _missingValue)
        : genofile(_genofile), idsfile(_idsfile), MissingValue(_missingValue), 
          num_samples(0), num_markers(0), use_subset(!_idsfile.empty())
    {
        process_large_genotypes(genofile, idsfile, MissingValue);
    }

    // Constructor para compatibilidade com código anterior
    Functions(const std::string &_genofile, int _missingValue)
        : genofile(_genofile), idsfile(""), MissingValue(_missingValue), 
          num_samples(0), num_markers(0), use_subset(false)
    {
        process_large_genotypes(genofile, "", MissingValue);
    }

    // Getter methods
    size_t get_num_samples() const { return num_samples; }
    size_t get_num_markers() const { return num_markers; }
    const std::vector<std::string> &get_ids() const { return ids_list; }
    const std::vector<std::vector<double>> &get_snps() const { return snps; }

    void make_add_matrix();
    void make_dom_matrix();

    // Matrix inversion methods
    std::vector<std::vector<double>> cholesky_inverse(const std::vector<std::vector<double>> &matrix);

    std::vector<std::vector<double>> conjugate_gradient_inverse(
        const std::vector<std::vector<double>> &A,
        double tol = 1e-8,
        int max_iter = 1000);

    // Métodos atualizados com opção iterativa
    std::vector<std::vector<double>> get_inverse_additive_matrix(
        bool use_iterative = true,
        double lambda = 0.01,
        double tol = 1e-8);

    std::vector<std::vector<double>> get_inverse_dominance_matrix(
        bool use_iterative = true,
        double lambda = 0.01,
        double tol = 1e-8);

    // Adicione estas declarações na seção pública da classe Functions
    std::vector<std::vector<double>> get_inverse_additive_additive_matrix(
        bool use_iterative = true,
        double lambda = 0.01,
        double tol = 1e-8);

    std::vector<std::vector<double>> get_inverse_dominance_dominance_matrix(
        bool use_iterative = true,
        double lambda = 0.01,
        double tol = 1e-8);

    std::vector<std::vector<double>> get_inverse_additive_dominance_matrix(
        bool use_iterative = true,
        double lambda = 0.01,
        double tol = 1e-8);

    // std::vector<std::vector<double>> get_inverse_identity_matrix();
    void save_vectorized(const std::vector<std::vector<double>> &matrix,
                         const std::string &filename,
                         bool one_based_indexing = true);

    template <typename T>
    void print_matrix(const std::vector<std::vector<T>> &matrix,
                      const std::string &title = "Matrix",
                      int max_size = 5)
    {
        std::cout << "\n"
                  << title << " (" << matrix.size() << "x" << matrix.size() << "):\n";
        std::cout << std::string(50, '-') << "\n";

        if (matrix.empty())
        {
            std::cout << "Matriz vazia!\n";
            return;
        }

        int rows_to_show = std::min(max_size, (int)matrix.size());
        int cols_to_show = std::min(max_size, (int)matrix[0].size());

        // Cabeçalho das colunas
        std::cout << std::setw(8) << " ";
        for (int j = 0; j < cols_to_show; ++j)
        {
            std::cout << std::setw(12) << "[" + std::to_string(j) + "]";
        }
        if (cols_to_show < (int)matrix[0].size())
        {
            std::cout << std::setw(8) << "...";
        }
        std::cout << "\n";

        // Linhas da matriz
        for (int i = 0; i < rows_to_show; ++i)
        {
            std::cout << std::setw(6) << "[" + std::to_string(i) + "]";

            for (int j = 0; j < cols_to_show; ++j)
            {
                std::cout << std::setw(12) << std::fixed << std::setprecision(6)
                          << matrix[i][j];
            }

            if (cols_to_show < (int)matrix[i].size())
            {
                std::cout << std::setw(8) << "...";
            }
            std::cout << "\n";
        }

        if (rows_to_show < (int)matrix.size())
        {
            std::cout << std::setw(6) << "...";
            for (int j = 0; j < cols_to_show; ++j)
            {
                std::cout << std::setw(12) << "...";
            }
            std::cout << "\n";
        }

        std::cout << std::string(50, '-') << "\n";
    }
};

#endif // MA_FUNCTIONS_H
