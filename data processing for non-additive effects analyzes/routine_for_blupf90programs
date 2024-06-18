import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings("ignore")

def leitura_dados(geno_file, data_file):
    # Ler os arquivos CSV
    geno = pd.read_csv(geno_file, delimiter=",", header=0)
    data = pd.read_csv(data_file, sep=" ", header=0)

    # Classificar os DataFrames pelo ID
    geno.sort_values(by="ID", inplace=True)
    data.sort_values(by="ID", inplace=True)

    return geno, data

def merge_and_format_files(geno, data, output_file=None):
    # Classificar os DataFrames pelo ID
    geno.sort_values(by="ID", inplace=True)
    data.sort_values(by="ID", inplace=True)

    # Encontrar IDs comuns
    ids_comuns = pd.merge(data[["ID"]], geno[["ID"]], how="inner", on="ID")

    # Manter apenas os IDs comuns em 'data'
    data = data[data["ID"].isin(ids_comuns["ID"])]

    # Adicionar nova coluna de ID numerado sequencialmente
    data["newID"] = range(1, len(data) + 1)
    geno["ID"] = range(1, len(data) + 1)

    data.to_csv("renum_data.txt", sep=" ", index=False)

    if output_file is not None and output_file.strip() != "":
        # Fazer uma cópia do DataFrame geno
        geno_copy = geno.copy()
        
        # Converter IDs para strings e encontrar o comprimento máximo do ID
        geno_copy["ID"] = geno_copy["ID"].astype(str)
        max_id_length = max(geno_copy["ID"].apply(len))
        
        # Converter DataFrame para o formato desejado
        geno_formatado = (
            geno_copy.set_index("ID").applymap(str).apply(lambda x: "".join(x), axis=1)
        )
        # Salvar os dados no formato desejado
        with open(output_file + ".txt", "w") as f:
            for index, row in geno_formatado.items():
                # Formatar o ID com espaço fixo
                formatted_id = index.ljust(max_id_length + 1)
                
                # Escrever no arquivo no formato fixo
                f.write(f"{formatted_id}{row}\n")


    if output_file is None or output_file.strip() == "":
        return geno
    else:
        return geno

def calcular_inversa(matrix):
    U, S, Vt = np.linalg.svd(matrix)
    S_inv = np.where(S < 1e-10, 0, 1 / S)
    matrix_inv = np.dot(Vt.T * S_inv, U.T)
    return matrix_inv


def matriz_aditiva_fast(df):
    ids = df.iloc[:, 0]
    # Substituindo os valores 5 por NaN
    df_values = df.iloc[:, 1:].replace(5, np.nan).values

    n = df_values.shape[0]
    p1 = np.nansum(df_values, axis=0) / (2 * n)  # Usando np.nansum para lidar com NaN
    q1 = 1 - p1
    Z = np.where(
        df_values == 2, 2 - 2 * p1, np.where(df_values == 1, 1 - 2 * p1, -2 * p1)
    )
    denominador = 2 * np.nansum(p1 * q1)  # Usando np.nansum para lidar com NaN

    GMATRIX = np.dot(Z, Z.T) / denominador

    # Calculando a pseudo-inversa
    GMATRIX_inv = calcular_inversa(GMATRIX)

    # Criando DataFrame com rótulos nas colunas e linhas
    GMATRIX_df = pd.DataFrame(GMATRIX_inv, index=ids, columns=ids)

    return GMATRIX_df


def matriz_dominancia_fast(df):

    ids = df.iloc[:, 0]
    df_values = df.iloc[:, 1:].replace(5, np.nan).values
    n = df_values.shape[0]
    p1 = np.nansum(df_values, axis=0) / (2 * n)
    q1 = 1 - p1
    W = np.where(
        df_values == 2,
        -2 * (q1**2),
        np.where(df_values == 1, 2 * (p1 * q1), -2 * (p1**2)),
    )
    denominador = np.nansum((2 * p1 * q1) ** 2)
    DMATRIX = np.dot(W, W.T) / denominador
    # Invertendo a matriz GMATRIX
    DMATRIX_invertida = calcular_inversa(DMATRIX)

    # Criando DataFrame com rótulos nas colunas e linhas
    DMATRIX_df = pd.DataFrame(DMATRIX_invertida, index=ids, columns=ids)

    return DMATRIX_df


def matriz_aditiva_aditiva_fast(df):
    ids = df.iloc[:, 0]
    # Substituindo os valores 5 por NaN
    df_values = df.iloc[:, 1:].replace(5, np.nan).values

    n = df_values.shape[0]
    p1 = np.nansum(df_values, axis=0) / (2 * n)  # Usando np.nansum para lidar com NaN
    q1 = 1 - p1
    Z = np.where(
        df_values == 2, 2 - 2 * p1, np.where(df_values == 1, 1 - 2 * p1, -2 * p1)
    )
    denominador = 2 * np.nansum(p1 * q1)  # Usando np.nansum para lidar com NaN

    GMATRIX = np.dot(Z, Z.T) / denominador
    # Normalização
    num = GMATRIX * GMATRIX
    denominador2 = np.trace(num) / n
    G_aditiva_aditiva_fast = num / denominador2

    GAAMATRIX_invertida = calcular_inversa(G_aditiva_aditiva_fast)

    # Criando DataFrame com rótulos nas colunas e linhas
    GAAMATRIX_df = pd.DataFrame(GAAMATRIX_invertida, index=ids, columns=ids)

    return GAAMATRIX_df


def matriz_dominancia_dominancia_fast(df):
    ids = df.iloc[:, 0]
    df_values = df.iloc[:, 1:].replace(5, np.nan).values
    n = df_values.shape[0]
    p1 = np.nansum(df_values, axis=0) / (2 * n)
    q1 = 1 - p1
    W = np.where(
        df_values == 2,
        -2 * (q1**2),
        np.where(df_values == 1, 2 * (p1 * q1), -2 * (p1**2)),
    )
    denominador = np.nansum((2 * p1 * q1) ** 2)
    DMATRIX = np.dot(W, W.T) / denominador
    # Normalização
    num = DMATRIX * DMATRIX
    denominador2 = np.trace(num) / n
    D_dominancia_dominancia = num / denominador2

    DDMATRIX_invertida = calcular_inversa(D_dominancia_dominancia)

    # Criando DataFrame com rótulos nas colunas e linhas
    DDMATRIX_df = pd.DataFrame(DDMATRIX_invertida, index=ids, columns=ids)

    return DDMATRIX_df


def matriz_aditiva_dominancia_fast(df):
    ids = df.iloc[:, 0]
    # Substituindo os valores 5 por NaN
    df_values = df.iloc[:, 1:].replace(5, np.nan).values

    n = df_values.shape[0]
    p1 = np.nansum(df_values, axis=0) / (2 * n)  # Usando np.nansum para lidar com NaN
    q1 = 1 - p1
    Z = np.where(
        df_values == 2, 2 - 2 * p1, np.where(df_values == 1, 1 - 2 * p1, -2 * p1)
    )
    denominador = 2 * np.nansum(p1 * q1)  # Usando np.nansum para lidar com NaN

    GMATRIX = np.dot(Z, Z.T) / denominador

    W = np.where(
        df_values == 2,
        -2 * (q1**2),
        np.where(df_values == 1, 2 * (p1 * q1), -2 * (p1**2)),
    )
    denominador = np.nansum((2 * p1 * q1) ** 2)
    DMATRIX = np.dot(W, W.T) / denominador

    # Normalização
    num = GMATRIX * DMATRIX
    denominador2 = np.trace(num) / n
    D_aditiva_dominancia = num / denominador2

    ADMATRIX_invertida = calcular_inversa(D_aditiva_dominancia)

    # Criando DataFrame com rótulos nas colunas e linhas
    ADMATRIX_df = pd.DataFrame(ADMATRIX_invertida, index=ids, columns=ids)

    return ADMATRIX_df


def salvar_vetorizado(matriz, nome_arquivo):
    # Obtém os índices das linhas e colunas
    ids_linhas = matriz.index
    ids_colunas = matriz.columns

    # Inicializa listas para armazenar os dados vetorizados
    vetorizado = []

    # Percorre todas as células da matriz
    for linha in ids_linhas:
        for coluna in ids_colunas:
            valor = matriz.loc[linha, coluna]  # Obtém o valor da célula
            if valor != 0 and coluna >= linha:  # Considera apenas valores diferentes de zero e coluna maior ou igual à linha
                vetorizado.append([linha, coluna, valor])  # Adiciona à lista vetorizada

    # Ordena os dados conforme especificado
    vetorizado.sort(key=lambda x: (x[0], x[1]))  # Ordena primeiro pela linha e depois pela coluna

    # Converte a lista vetorizada em um DataFrame do Pandas
    vetorizado_df = pd.DataFrame(vetorizado, columns=['Linha', 'Coluna', 'Valor'])

    # Salva o DataFrame em um arquivo de texto
    vetorizado_df.to_csv(nome_arquivo, sep=' ', index=False, float_format='%.16f', header=False)

def call_rate(df, call_rate_min_snps, call_rate_min_samples, missing_value):
    # Dicionário para armazenar as proporções de ocorrência do valor 5 em cada coluna
    missing_snps = {}
    # Iterar sobre todas as colunas do DataFrame, começando da segunda coluna (índice 1)
    for coluna in df.columns[1:]:
        # Contar o número de ocorrências do valor 5 na coluna atual
        ocorrencias_missing = np.sum(df[coluna] == missing_value)

        # Calcular a proporção de ocorrências do valor 5 em relação ao tamanho da coluna
        n_missing = ocorrencias_missing / len(df[coluna])

        # Armazenar a proporção no dicionário
        missing_snps[coluna] = n_missing

    # Identificar SNPs para serem removidos (onde a proporção de missing_value é maior que call_rate_max)
    snps_call_rate = [
        coluna
        for coluna, call_rate in missing_snps.items()
        if 1 - call_rate < call_rate_min_snps
    ]
    # Remover os SNPs identificados
    df_limpo_snps = df.drop(columns=snps_call_rate)

    # Mensagem de impressão informando ao usuário quais SNPs foram removidos
    num_snps_removidos = len(snps_call_rate)

    # Mensagem de impressão informando ao usuário o número de SNPs removidos
    if num_snps_removidos > 0:
        print(f"{num_snps_removidos} SNPs were removed due to low call rate:")
        for coluna in snps_call_rate:
            print(coluna)
    else:
        print("No SNPs were removed due to low call rate.")
    # Lista para armazenar os índices e IDs das amostras com baixo call rate
    linhas_remover = []

    # Iterar sobre todas as linhas do DataFrame, começando da segunda coluna (índice 1)
    for indice, linha in df_limpo_snps.iloc[:, 1:].iterrows():
        # Contar o número de ocorrências do valor 5 na linha atual
        ocorrencias_missing = np.nansum(linha == missing_value)

        # Calcular a proporção de ocorrências do valor 5 em relação ao tamanho da linha
        n_missing = ocorrencias_missing / len(linha)

        # Verificar se a proporção é menor que call_rate_min
        if n_missing > 1 - call_rate_min_samples:
            linhas_remover.append(indice)

    # Remover as linhas identificadas
    df_limpo = df_limpo_snps.drop(index=linhas_remover)
    print()
    print()
    # Mensagem de impressão informando ao usuário quantas linhas foram removidas e quais IDs foram removidos
    if linhas_remover:
        print(
            f"{len(linhas_remover)} Samples were removed due to low call rate:"
        )
        for indice in linhas_remover:
            print(f"ID: {df.loc[indice, 'ID']}")
    else:
        print("No samples were removed due to low call rate.")
    return df_limpo


def maf(df, min_maf=0.05):
    # Calcular as frequências alélicas
    df_values = df.values[:, 1:]
    n = df_values.shape[0]
    p1 = np.sum(df_values, axis=0) / (2 * n)
    # Identificar colunas para serem removidas (onde a frequência alélica é menor que min_maf)
    cols_to_remove = [col for col, freq in enumerate(p1) if freq < min_maf]

    # Remover as colunas identificadas
    df_limpo = df.drop(columns=df.columns[cols_to_remove])

    # Mensagem de impressão informando ao usuário quais colunas foram removidas
    if cols_to_remove:
        print(
            "The following columns have been removed due to low Minimum Allele Frequency (MAF):"
        )
        for col in cols_to_remove:
            print(df.columns[col])
    else:
        print(
            "No columns were removed due to low Minimum Allele Frequency (MAF)."
        )

    return df_limpo


def monomorphic_snp(df, missing_value=5):
    # Verificar se todos os valores em cada coluna são iguais, exceto os valores de missing
    monomorphic_cols = []
    for coluna in df.columns[1:]:
        unique_values = df[coluna][df[coluna] != missing_value].unique()
        if len(unique_values) == 1:
            monomorphic_cols.append(coluna)

    # Remover as colunas identificadas
    df_limpo = df.drop(columns=monomorphic_cols)

    # Mensagem de impressão informando ao usuário quais colunas foram removidas
    if monomorphic_cols:
        print("The following columns were removed due to monomorphic markers:")
        for col in monomorphic_cols:
            print(col)
    else:
        print("No columns were removed due to monomorphic markers.")

    return df_limpo


def hwe_test(df, limiar=0.15):
    # Calcular as frequências alélicas
    df_values = df.values[:, 1:]
    n = df_values.shape[0]
    p1 = np.sum(df_values, axis=0) / (2 * n)
    q1 = 1 - p1
    freq_het = np.sum(df_values == 1, axis=0) / n

    # Calcular a frequência esperada de heterozigotos e o valor do teste de HWE para cada alelo em cada SNP
    cols_to_remove = []
    num_removed = 0  # Contador de marcadores removidos

    for i, coluna in enumerate(df.columns[1:]):
        freq_het_exp = 2 * p1[i] * q1[i]
        hwe_value = freq_het[i] - freq_het_exp
        # Verificar se o valor do teste de HWE é maior que 0.15 e marcar a coluna para remoção
        if hwe_value > limiar:
            cols_to_remove.append(coluna)
            num_removed += 1

    # Remover as colunas identificadas
    df_limpo = df.drop(columns=cols_to_remove)
    # Mensagem de impressão informando ao usuário quais colunas foram removidas
    if cols_to_remove:
        print(
            f"{num_removed} markers were removed due to the Hardy-Weinberg test:"
        )
        for col in cols_to_remove:
            print(col)
    else:
        print("No columns were removed due to the Hardy-Weinberg test.")

    return df_limpo


def detect_duplicates(df):
    # Verificar duplicatas
    duplicates = df[df.duplicated()]

    # Se houver duplicatas, imprimir e remover
    if not duplicates.empty:
        print("The following duplicate lines were detected and removed:")
        print(duplicates)
        df = df.drop_duplicates()
    else:
        print("No duplicate lines were found.")

    return df


def quality_control(
    df,
    call_rate_min_snps=0.90,
    call_rate_min_samples=0.90,
    min_maf=0.05,
    limiar_hwe=0.15,
    missing_value=5,
):
    # Testar para Call Rate
    print("Testing for Call Rate...")
    df = call_rate(df, call_rate_min_snps, call_rate_min_samples, missing_value)
    print()

    # Testar para MAF
    print("Testing for MAF...")
    df = maf(df, min_maf)
    print()

    # Testar para Monomórfico
    print("Testing for monomorphic...")
    df = monomorphic_snp(df, missing_value)
    print()

    # Testar para Hardy-Weinberg Equilíbrio
    print("Testing for Hardy-Weinberg test...")
    df = hwe_test(df, limiar_hwe)
    print()

    # Detectar e remover duplicatas
    print("Detecting and removing duplicates...")
    df = detect_duplicates(df)
    print()


    return df

# Chamando a função de controle de qualidade


def gerar_matriz(df, tipo_matriz, nome_matriz):
    if tipo_matriz == "Additive":
        matriz = matriz_aditiva_fast(df)
    elif tipo_matriz == "Dominance":
        matriz = matriz_dominancia_fast(df)
    elif tipo_matriz == "Additive_Additive":
        matriz = matriz_aditiva_aditiva_fast(df)
    elif tipo_matriz == "Dominance_Dominance":
        matriz = matriz_dominancia_dominancia_fast(df)
    elif tipo_matriz == "Additive_Dominance":
        matriz = matriz_aditiva_dominancia_fast(df)
    else:
        raise ValueError(
            "Invalid array type. Choose between 'additive', 'dominance', 'additive by additive', 'dominance by dominance' or 'additive by dominance''."
        )
    print()
    print("Matrix ready and inverted, wait to vectorize...")
    nome_arquivo = nome_matriz + ".txt"
    salvar_vetorizado(matriz, nome_arquivo)
    print()
    print("Vectorized matrix  =D  ")


def main():
    while True:
        # Solicitação dos nomes dos arquivos de entrada e saída ao usuário
        print("It is mandatory to have a column identified as 'ID' in this way, geno must be a tabular form also")
        geno_file = input(r"Enter the path of the CSV file containing the geno file: ")
        data_file = input(r"Enter the path of the CSV file containing the data file: ")
        output_file = input(
            "Enter the name of the formatted output file (or leave it blank to not save), save as blupf90 format: "
        )

        # Leitura dos dados
        geno_df, data_df = leitura_dados(geno_file, data_file)

        # Exibição do cabeçalho dos dados
        print("The first rows of your geno file")
        print(geno_df.head(4))


        # Solicitação de entrada do usuário para controle de qualidade
        controle_qualidade = input(
            "Do you want to perform data quality control? (y/n): "
        )

        # Verificação da entrada do usuário e execução do controle de qualidade
        if controle_qualidade.lower() == "y":
            df_limpo = quality_control(geno_df)
        elif controle_qualidade.lower() == "n":
            df_limpo = geno_df
        else:
            print("Invalid Input. Please provide 'y' for yes or 'n' for no.")
            continue
        # Merge e formatação dos dados
        geno_df = merge_and_format_files(df_limpo, data_df, output_file)

        while True:
            # Solicitação do tipo de matriz e nome do arquivo
            print("Available matrix types:")
            print("1 - Additive")
            print("2 - Dominance")
            print("3 - Additive by Additive ")
            print("4 - Dominance by Dominance")
            print("5 - Additive Dominance")
            print("6 - Choose another dataset")
            print("7 - Stop the program")
            print()
            print()
            tipo_matriz_numero = int(
                input("Choose the number corresponding to the desired matrix type: ")
            )
            print()
            if tipo_matriz_numero == 7:
                print("Program closed.")
                return
            elif tipo_matriz_numero == 6:
                break

            nome_matriz = input(
                "Enter the name for the matrix to be generated (without extension): "
            )

            # Mapear o número escolhido para o tipo de matriz correspondente
            if tipo_matriz_numero == 1:
                tipo_matriz = "Additive"
            elif tipo_matriz_numero == 2:
                tipo_matriz = "Dominance"
            elif tipo_matriz_numero == 3:
                tipo_matriz = "Additive_Additive"
            elif tipo_matriz_numero == 4:
                tipo_matriz = "Dominance_Dominance"
            elif tipo_matriz_numero == 5:
                tipo_matriz = "Additive_Dominance"
            else:
                print("Opção inválida. Por favor, escolha uma opção válida.")
                continue

            # Geração da matriz
            gerar_matriz(geno_df, tipo_matriz, nome_matriz)

        continuar = input("Do you want to generate another matrix? (y/n): ")
        if continuar.lower() != "y":
            print("Program closed")
            break


# Execução da função main
if __name__ == "__main__":
    main()
