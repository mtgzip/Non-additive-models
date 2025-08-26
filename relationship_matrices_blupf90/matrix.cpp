#include "matrix.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

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

void Functions::load_subset_ids(const std::string &ids_file)
{
    if (ids_file.empty())
    {
        use_subset = false;
        return;
    }

    std::ifstream file(ids_file);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open IDs file: " + ids_file);
    }

    // Limpa estruturas anteriores
    ordered_subset_ids.clear();
    subset_ids.clear();
    
    std::string line;
    size_t count = 0;

    std::cout << " " << std::endl;

    std::cout << "Loading subset IDs from: " << ids_file << std::endl;


    while (std::getline(file, line))
    {
        // Remove espaços em branco no início e fim
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (!line.empty())
        {
            ordered_subset_ids.push_back(line);  // Mantém a ordem
            subset_ids.insert(line);             // Para busca rápida
            count++;
        }
    }

    file.close();
    use_subset = true;

    std::cout << "Loaded " << count << " IDs for subsetting in specified order" << std::endl;

    if (ordered_subset_ids.empty())
    {
        std::cout << "Warning: No valid IDs found in subset file. Processing all individuals." << std::endl;
        use_subset = false;
    }
}

void Functions::process_large_genotypes(const std::string &file_path, const std::string &ids_file, int missingValue)
{
    // Carrega IDs do subset se fornecido
    load_subset_ids(ids_file);

    std::ifstream file(file_path);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open file: " + file_path);
    }

    // Primeiro passo: carregar todos os genótipos em um mapa temporário
    std::cout << "Loading genotype data into temporary storage..." << std::endl;
        std::cout << " " << std::endl;

    std::unordered_map<std::string, std::vector<double>> genotype_map;
    size_t n_markers = 0;
    bool first_line = true;
    
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string id, genotype_str;

        if (iss >> id >> genotype_str)
        {
            // Se estamos usando subset, verifica se o ID está na lista
            if (use_subset && subset_ids.find(id) == subset_ids.end())
            {
                continue; // Pula este indivíduo
            }

            if (first_line)
            {
                n_markers = genotype_str.size();
                first_line = false;
            }

            // Processa genótipos
            std::vector<double> individual_genotypes;
            individual_genotypes.reserve(n_markers);

            for (size_t j = 0; j < std::min(genotype_str.size(), n_markers); j++)
            {
                char c = genotype_str[j];
                if (c >= '0' && c <= '9')
                {
                    individual_genotypes.push_back(c - '0');
                }
                else
                {
                    individual_genotypes.push_back(missingValue);
                }
            }

            // Preenche com missing values se necessário
            while (individual_genotypes.size() < n_markers)
            {
                individual_genotypes.push_back(missingValue);
            }

            genotype_map[id] = std::move(individual_genotypes);
        }
    }
    file.close();

    // Segundo passo: organizar na ordem correta
    std::cout << "Organizing data in specified order..." << std::endl;
        std::cout << " " << std::endl;

    ids_list.clear();
    snps.clear();
    
    if (use_subset)
    {
        // Usar a ordem dos IDs do arquivo de subset
        ids_list.reserve(ordered_subset_ids.size());
        snps.reserve(ordered_subset_ids.size());
        
        size_t found_count = 0;
        for (const std::string& id : ordered_subset_ids)
        {
            auto it = genotype_map.find(id);
            if (it != genotype_map.end())
            {
                ids_list.push_back(id);
                snps.push_back(std::move(it->second));
                found_count++;
            }
            else
            {
                std::cout << "Warning: ID '" << id << "' from subset file not found in genotype file" << std::endl;
            }
        }
        
        std::cout << "Found " << found_count << " out of " << ordered_subset_ids.size() << " requested IDs" << std::endl;
    }
    else
    {
        // Se não há subset, usar todos na ordem que aparecem no arquivo
        std::ifstream file2(file_path);
        std::string line2;
        while (std::getline(file2, line2))
        {
            if (line2.empty()) continue;
            
            std::istringstream iss(line2);
            std::string id, genotype_str;
            
            if (iss >> id >> genotype_str)
            {
                auto it = genotype_map.find(id);
                if (it != genotype_map.end())
                {
                    ids_list.push_back(id);
                    snps.push_back(std::move(it->second));
                }
            }
        }
        file2.close();
    }

    // Atualiza dimensões
    num_samples = ids_list.size();
    num_markers = n_markers;

    std::cout << "File processing completed successfully!" << std::endl;
    std::cout << "Final dataset: " << num_samples << " samples with " << num_markers << " markers each." << std::endl;
    
    // Mostra os primeiros IDs na ordem final
    std::cout << "\nFinal ID order (first 10):" << std::endl;
    for (size_t i = 0; i < std::min(ids_list.size(), size_t(10)); i++)
    {
        std::cout << "  [" << i << "] " << ids_list[i] << std::endl;
    }
}

template <typename T>
void compute_symmetric_product_and_normalize(
    const std::vector<std::vector<T>> &Z,
    std::vector<std::vector<T>> &result,
    T normalizer)
{

    size_t num_rows = Z.size();
    size_t num_cols = Z[0].size();

    // Redimensiona a matriz de resultado, se necessário
    result.resize(num_rows, std::vector<T>(num_rows, T()));

    // Calcula Z*Z^T aproveitando a simetria
    for (size_t i = 0; i < num_rows; ++i)
    {
        for (size_t j = i; j < num_rows; ++j)
        {
            T sum = 0;
            for (size_t k = 0; k < num_cols; ++k)
            {
                sum += Z[i][k] * Z[j][k];
            }
            // Normaliza diretamente
            result[i][j] = sum / normalizer;
            if (i != j)
            {
                result[j][i] = result[i][j]; // Aproveita a simetria
            }
        }
    }
}

void Functions::make_add_matrix()
{
    // 1. Calcular frequências alélicas (p) e denominador
    std::vector<double> p(num_markers, 0);
    std::vector<int> valid_counts(num_markers, 0);
    double denominator = 0.0;

    for (size_t col = 0; col < num_markers; ++col)
    {
        // Calcular frequência alélica para SNPs válidos
        for (size_t row = 0; row < num_samples; ++row)
        {
            int val = snps[row][col];
            if (val != MissingValue)
            {
                p[col] += val;
                valid_counts[col]++;
            }
        }

        // Calcular frequência e contribuição para o denominador
        if (valid_counts[col] > 0)
        {
            p[col] /= (2.0 * valid_counts[col]);
            denominator += 2.0 * p[col] * (1.0 - p[col]);
        }
        else
        {
            p[col] = -1.0; // SNP inválido
        }
    }

    // 2. Construir matriz Z centrada
    std::vector<std::vector<double>> Z(num_samples, std::vector<double>());
    size_t valid_snps = 0;

    for (size_t row = 0; row < num_samples; ++row)
    {
        Z[row].reserve(num_markers);
        for (size_t col = 0; col < num_markers; ++col)
        {
            // Incluir apenas SNPs válidos e não excluídos

            Z[row].push_back(snps[row][col] - 2.0 * p[col]);
            valid_snps++;
        }
    }

    // 3. Computar matriz de relacionamento aditivo
    additive_matrix.clear();
    if (valid_snps > 0)
    {
        compute_symmetric_product_and_normalize(Z, additive_matrix, denominator);
    }
    else
    {
        additive_matrix.resize(num_samples, std::vector<double>(num_samples, 0.0));
    }

    // // Debug (opcional)

    print_matrix(additive_matrix, "G (Additive)");
}

void Functions::make_dom_matrix()
{
    // 1. Calcular frequências alélicas (p) e denominador
    std::vector<double> p(num_markers, 0);
    std::vector<int> valid_counts(num_markers, 0);
    double denominator = 0.0;

    for (size_t col = 0; col < num_markers; ++col)
    {

        // Calcular frequência alélica para SNPs válidos
        for (size_t row = 0; row < num_samples; ++row)
        {
            int val = snps[row][col];
            if (val != MissingValue)
            {
                p[col] += val;
                valid_counts[col]++;
            }
        }

        // Calcular frequência e contribuição para o denominador
        if (valid_counts[col] > 0)
        {
            p[col] /= (2.0 * valid_counts[col]);
            double het = 2.0 * p[col] * (1.0 - p[col]);
            denominator += het * het; // Diferente da matriz aditiva!
        }
        else
        {
            p[col] = -1.0; // SNP inválido
        }
    }

    // 2. Construir matriz W de dominância
    std::vector<std::vector<double>> W(num_samples, std::vector<double>());
    size_t valid_snps = 0;

    for (size_t row = 0; row < num_samples; ++row)
    {
        W[row].reserve(num_markers);
        for (size_t col = 0; col < num_markers; ++col)
        {
            // Incluir apenas SNPs válidos e não excluídos

            double p_val = p[col];
            double q = 1.0 - p_val;
            int geno = snps[row][col];

            // Codificação de dominância
            if (geno == 2)
            {
                W[row].push_back(-2.0 * q * q);
            }
            else if (geno == 1)
            {
                W[row].push_back(2.0 * p_val * q);
            }
            else if (geno == 0)
            {
                W[row].push_back(-2.0 * p_val * p_val);
            }
            else
            {
                W[row].push_back(0.0); // Missing
            }
            valid_snps++;
        }
    }

    // 3. Computar matriz de relacionamento de dominância
    dominance_matrix.clear();
    if (valid_snps > 0)
    {
        compute_symmetric_product_and_normalize(W, dominance_matrix, denominator);
    }
    else
    {
        dominance_matrix.resize(num_samples, std::vector<double>(num_samples, 0.0));
    }

    // Debug (opcional)

    print_matrix(dominance_matrix, "D (Dominance)");
}


std::vector<std::vector<double>> Functions::cholesky_decomposition(const std::vector<std::vector<double>> &matrix)
{
    size_t n = matrix.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j <= i; j++)
        {
            double sum = 0.0;

            if (j == i)
            {
                for (size_t k = 0; k < j; k++)
                {
                    sum += L[j][k] * L[j][k];
                }
                double diag = matrix[j][j] - sum;
                if (diag <= 0)
                {
                    // Adiciona tolerância numérica para diagonais não-positivas
                    diag = 1e-8;
                }
                L[j][j] = std::sqrt(diag);
            }
            else
            {
                for (size_t k = 0; k < j; k++)
                {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (matrix[i][j] - sum) / L[j][j];
            }
        }
    }

    return L;
}

// Forward Substitution (L * Y = B)
std::vector<std::vector<double>> Functions::forward_substitution(const std::vector<std::vector<double>> &L, const std::vector<std::vector<double>> &B)
{
    size_t n = L.size();
    size_t m = B[0].size();
    std::vector<std::vector<double>> Y(n, std::vector<double>(m, 0.0));

    for (size_t j = 0; j < m; j++)
    {
        for (size_t i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (size_t k = 0; k < i; k++)
            {
                sum += L[i][k] * Y[k][j];
            }
            Y[i][j] = (B[i][j] - sum) / L[i][i];
        }
    }

    return Y;
}

// Backward Substitution (L^T * X = Y)
std::vector<std::vector<double>> Functions::backward_substitution(const std::vector<std::vector<double>> &U, const std::vector<std::vector<double>> &Y)
{
    size_t n = U.size();
    size_t m = Y[0].size();
    std::vector<std::vector<double>> X(n, std::vector<double>(m, 0.0));

    for (size_t j = 0; j < m; j++)
    {
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0.0;
            for (size_t k = i + 1; k < n; k++)
            {
                sum += U[k][i] * X[k][j]; // Note: U[k][i] porque U = L^T
            }
            X[i][j] = (Y[i][j] - sum) / U[i][i];
        }
    }

    return X;
}

// Matrix Multiplication
std::vector<std::vector<double>> Functions::matrix_multiply(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B)
{
    size_t n = A.size();
    size_t m = B[0].size();
    size_t p = B.size();

    std::vector<std::vector<double>> C(n, std::vector<double>(m, 0.0));

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            for (size_t k = 0; k < p; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

// Matrix Transpose
std::vector<std::vector<double>> Functions::matrix_transpose(const std::vector<std::vector<double>> &matrix)
{
    size_t n = matrix.size();
    size_t m = matrix[0].size();

    std::vector<std::vector<double>> transposed(m, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

// Cholesky Inverse
std::vector<std::vector<double>> Functions::cholesky_inverse(const std::vector<std::vector<double>> &matrix)
{
    // Verifica se a matriz é quadrada
    if (matrix.size() != matrix[0].size())
    {
        throw std::invalid_argument("Matrix must be square for inversion");
    }

    size_t n = matrix.size();

    // 1. Decomposição de Cholesky: A = L * L^T
    std::vector<std::vector<double>> L = cholesky_decomposition(matrix);

    // 2. Cria matriz identidade
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; i++)
    {
        I[i][i] = 1.0;
    }

    // 3. Resolve L * Y = I usando substituição direta
    std::vector<std::vector<double>> Y = forward_substitution(L, I);

    // 4. Resolve L^T * X = Y usando substituição reversa
    std::vector<std::vector<double>> X = backward_substitution(L, Y);

    return X;
}

std::vector<std::vector<double>> Functions::conjugate_gradient_inverse(
    const std::vector<std::vector<double>> &A,
    double tol,
    int max_iter)
{
    size_t n = A.size();
    std::vector<std::vector<double>> invA(n, std::vector<double>(n, 0.0));

    // Para cada coluna da matriz identidade
    for (size_t j = 0; j < n; j++)
    {
        // Vetor b = j-ésima coluna da identidade
        std::vector<double> b(n, 0.0);
        b[j] = 1.0;

        // Vetor solução x
        std::vector<double> x(n, 0.0);
        std::vector<double> r = b; // resíduo inicial
        std::vector<double> p = r; // direção de busca

        double r_norm_old = 0.0;
        for (size_t i = 0; i < n; i++)
        {
            r_norm_old += r[i] * r[i];
        }

        int iter = 0;
        while (iter < max_iter)
        {
            // Calcula A*p
            std::vector<double> Ap(n, 0.0);
            for (size_t i = 0; i < n; i++)
            {
                for (size_t k = 0; k < n; k++)
                {
                    Ap[i] += A[i][k] * p[k];
                }
            }

            // alpha = (r·r) / (p·A·p)
            double pAp = 0.0;
            for (size_t i = 0; i < n; i++)
            {
                pAp += p[i] * Ap[i];
            }
            double alpha = r_norm_old / pAp;

            // Atualiza x e r
            for (size_t i = 0; i < n; i++)
            {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }

            // Verifica convergência
            double r_norm_new = 0.0;
            for (size_t i = 0; i < n; i++)
            {
                r_norm_new += r[i] * r[i];
            }

            if (std::sqrt(r_norm_new) < tol)
            {
                break;
            }

            // Atualiza direção de busca
            double beta = r_norm_new / r_norm_old;
            for (size_t i = 0; i < n; i++)
            {
                p[i] = r[i] + beta * p[i];
            }

            r_norm_old = r_norm_new;
            iter++;
        }

        // Armazena a coluna j da inversa
        for (size_t i = 0; i < n; i++)
        {
            invA[i][j] = x[i];
        }
    }

    return invA;
}
std::vector<std::vector<double>> Functions::get_inverse_additive_matrix(
    bool use_iterative, double lambda, double tol)
{
    if (additive_matrix.empty())
        make_add_matrix();

    std::vector<std::vector<double>> A_reg = additive_matrix;
    for (size_t i = 0; i < A_reg.size(); i++)
    {
        A_reg[i][i] += lambda;
    }

    if (use_iterative)
    {
        return conjugate_gradient_inverse(A_reg, tol);
    }
    else
    {
        return cholesky_inverse(A_reg);
    }
}

std::vector<std::vector<double>> Functions::get_inverse_dominance_matrix(
    bool use_iterative, double lambda, double tol)
{
    if (dominance_matrix.empty())
        make_dom_matrix();

    std::vector<std::vector<double>> D_reg = dominance_matrix;
    for (size_t i = 0; i < D_reg.size(); i++)
    {
        D_reg[i][i] += lambda;
    }

    if (use_iterative)
    {
        return conjugate_gradient_inverse(D_reg, tol);
    }
    else
    {
        return cholesky_inverse(D_reg);
    }
}

std::vector<std::vector<double>> Functions::get_inverse_additive_additive_matrix(
    bool use_iterative, double lambda, double tol)
{
    if (additive_matrix.empty())
    {
        make_add_matrix();
    }

    size_t n = additive_matrix.size();
    std::vector<std::vector<double>> GAA(n, std::vector<double>(n, 0.0));

    // Calcular G⊙G (element-wise square)
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            GAA[i][j] = additive_matrix[i][j] * additive_matrix[i][j];
        }
    }

    // Normalização: traço de GAA ao quadrado dividido por n
    double trace_sq = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        trace_sq += GAA[i][i];

        // trace_sq += GAA[i][i] * GAA[i][i];
    }
    double denominator = trace_sq / n;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            GAA[i][j] /= denominator;
        }
    }

    // Regularização
    for (size_t i = 0; i < n; i++)
    {
        GAA[i][i] += lambda;
    }

    if (use_iterative)
    {
        return conjugate_gradient_inverse(GAA, tol);
    }
    else
    {
        return cholesky_inverse(GAA);
    }
}

std::vector<std::vector<double>> Functions::get_inverse_dominance_dominance_matrix(
    bool use_iterative, double lambda, double tol)
{
    if (dominance_matrix.empty())
    {
        make_dom_matrix();
    }

    size_t n = dominance_matrix.size();
    std::vector<std::vector<double>> DDD(n, std::vector<double>(n, 0.0));

    // Calcular D⊙D (element-wise square)
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            DDD[i][j] = dominance_matrix[i][j] * dominance_matrix[i][j];
        }
    }

    // Normalização: traço de DDD ao quadrado dividido por n
    double trace_sq = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        trace_sq += DDD[i][i];
        // trace_sq += DDD[i][i] * DDD[i][i];
    }
    double denominator = trace_sq / n;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            DDD[i][j] /= denominator;
        }
    }

    // Regularização
    for (size_t i = 0; i < n; i++)
    {
        DDD[i][i] += lambda;
    }

    if (use_iterative)
    {
        return conjugate_gradient_inverse(DDD, tol);
    }
    else
    {
        return cholesky_inverse(DDD);
    }
}

std::vector<std::vector<double>> Functions::get_inverse_additive_dominance_matrix(
    bool use_iterative, double lambda, double tol)
{
    if (additive_matrix.empty())
    {
        make_add_matrix();
    }
    if (dominance_matrix.empty())
    {
        make_dom_matrix();
    }

    size_t n = additive_matrix.size();
    std::vector<std::vector<double>> ADD(n, std::vector<double>(n, 0.0));

    // Calcular G⊙D (element-wise product)
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            ADD[i][j] = additive_matrix[i][j] * dominance_matrix[i][j];
        }
    }

    // Normalização: traço de ADD ao quadrado dividido por n
    double trace_sq = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        trace_sq += ADD[i][i];
        // trace_sq += ADD[i][i] * ADD[i][i];
    }
    double denominator = trace_sq / n;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            ADD[i][j] /= denominator;
        }
    }

    // Regularização
    for (size_t i = 0; i < n; i++)
    {
        ADD[i][i] += lambda;
    }

    if (use_iterative)
    {
        return conjugate_gradient_inverse(ADD, tol);
    }
    else
    {
        return cholesky_inverse(ADD);
    }
}

void Functions::save_vectorized(const std::vector<std::vector<double>> &matrix,
                                const std::string &filename,
                                bool one_based_indexing)
{
    if (matrix.empty())
    {
        throw std::runtime_error("Matrix is empty");
    }

    size_t n = matrix.size();

    // Abre o arquivo imediatamente para streaming
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    outfile << std::scientific << std::setprecision(16);

    size_t count = 0;

    // Processa e salva diretamente, sem armazenar tudo na memória
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i; j < n; j++)
        {
            double value = matrix[i][j];
            if (value != 0.0)
            {
                size_t row = i;
                size_t col = j;

                if (one_based_indexing)
                {
                    row++;
                    col++;
                }

                outfile << row << " " << col << " " << value << "\n";
                count++;
            }
        }
    }

    outfile.close();
    std::cout << "Matrix saved to: " << filename
              << " (" << count << " non-zero elements)" << std::endl;
}
