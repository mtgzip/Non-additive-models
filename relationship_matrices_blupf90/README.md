# Genomic Relationship Matrices Calculator

A high-performance C++ tool for computing genomic relationship matrices and their inverses from SNP genotype data. Designed for quantitative genetics and genomic prediction applications.

## Features

- **Multiple Matrix Types**: Computes additive, dominance, and epistatic relationship matrices
- **Efficient Matrix Inversion**: Both direct (Cholesky) and iterative (Conjugate Gradient) methods
- **ID Subset Support**: Process specific individuals while maintaining desired order
- **Memory Efficient**: Handles large datasets with optimized memory usage
- **Flexible Output**: Sparse matrix format for efficient storage

## Supported Matrices

- **Additive (G)**: Standard genomic relationship matrix
- **Dominance (D)**: Dominance relationship matrix  
- **Additive×Additive (G#G)**: Additive epistatic interactions
- **Dominance×Dominance (D#D)**: Dominance epistatic interactions
- **Additive×Dominance (G#D)**: Additive-dominance epistatic interactions

## Installation

### Requirements
- C++11 or later
- GCC compiler

### Compilation
```bash
g++ -O3 -funroll-loops main.cpp matrix.cpp -o make_matrix
```

**Note**: This version is currently sequential. Future releases will include parallel implementations using OpenMP/threads and CUDA for GPU acceleration.

### Optional Optimizations
For better performance with large matrices:
```bash
g++ -O3 -funroll-loops -fopenmp main.cpp matrix.cpp -o make_matrix
```

## Usage

### Basic Usage
```bash
# Compute all matrices
./make_matrix genotype_file.txt

# Compute specific matrices
./make_matrix genotype_file.txt --matrices additive,dominance

# Use ID subset with custom missing value
./make_matrix genotype_file.txt --ids subset.txt --missing 9
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--ids <file>` | ID subset file (one ID per line) | None |
| `--missing <value>` | Missing genotype value | 5 |
| `--matrices <list>` | Comma-separated matrix types | all |
| `--help, -h` | Show help message | - |

### Matrix Types
- `additive`: Additive relationship matrix (G⁻¹)
- `dominance`: Dominance relationship matrix (D⁻¹)  
- `additive_additive`: Additive epistatic matrix (G#G)⁻¹
- `dominance_dominance`: Dominance epistatic matrix (D#D)⁻¹
- `additive_dominance`: Additive-dominance epistatic matrix (G#D)⁻¹
- `all`: Compute all available matrices

## Input Format

### Genotype File
Space-separated format with ID and genotype string:
```
ID1 012201210120...
ID2 120102012101...
ID3 201012101201...
```

**Genotype Encoding:**
- `0`: Homozygous reference (AA)
- `1`: Heterozygous (AB)  
- `2`: Homozygous alternate (BB)
- Other values: Treated as missing data

### ID Subset File (Optional)
One ID per line in desired output order:
```
ID1
ID3
ID5
```

**Important**: When using `--ids`, the final matrix order will match the order specified in the ID file, not the genotype file order.

## Output Files

The program generates sparse matrix files in the format:
```
row col value
1 1 1.245e+00
1 2 -3.421e-01
2 2 1.156e+00
...
```

### Output Files by Matrix Type
- `gi.txt`: Additive matrix inverse (G⁻¹)
- `di.txt`: Dominance matrix inverse (D⁻¹)
- `ggi.txt`: Additive×Additive matrix inverse
- `ddi.txt`: Dominance×Dominance matrix inverse  
- `gdi.txt`: Additive×Dominance matrix inverse

## Examples

### Example 1: Basic Usage
```bash
./make_matrix snp_data.txt --matrices additive
```
Computes only the additive relationship matrix inverse.

### Example 2: Subset Analysis
```bash
./make_matrix large_dataset.txt --ids selected_animals.txt --matrices all
```
Processes only animals listed in `selected_animals.txt` and computes all matrix types.

### Example 3: Custom Parameters
```bash
./make_matrix genotypes.txt --missing 9 --matrices additive,dominance
```
Uses `9` as missing value code and computes additive and dominance matrices only.

## Performance Considerations

### Memory Requirements
For a dataset with N individuals:
- **RAM needed**: ~8N² bytes (for double precision matrices)
- **Example**: 10,000 individuals ≈ 800 MB RAM

**Current Implementation**: Sequential processing
**Future Releases**: Will include parallel implementations using:
- **OpenMP/Threads**: Multi-core CPU parallelization
- **CUDA**: GPU acceleration for massive datasets

### Optimization Tips
1. **Large datasets**: The program automatically uses iterative methods for better performance
2. **Memory**: Ensure sufficient RAM for your dataset size
3. **CPU**: Multi-core processors will provide significant improvements in future parallel versions

## Algorithm Details

### Additive Matrix (G)
Based on VanRaden's (2008) method:
```
G = ZZ'/2Σ[pq)]
```
Where Z is the centered marker matrix and p is the allele frequency.

### Dominance Matrix (D)  
Following Vitezika et al. (2017) approach:
```
D = WW'/4Σ[p²q²]
```
Where W has elements equal -2q², 2pq, -2p² for 0,1,2 genotypes 

### Matrix Inversion
- **Cholesky decomposition**: For smaller, well-conditioned matrices
- **Conjugate Gradient**: For larger matrices or ill-conditioned cases
- **Regularization**: Adds small value (λ=0.01) to diagonal for numerical stability

## Citation

If you use this software in your research, please cite:

VanRaden, P.M. 2008. Efficient methods to compute genomic predictions. J. Dairy Sci. 91:4414-4423.

## License

This project is open source. Please check the repository for specific license terms.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Support

For questions or issues, please open an issue on the GitHub repository.
