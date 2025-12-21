# OMICS PROJECT - Biomarcadores em Cancro do Ovário

Este repositório contém a reanálise dos dados proteómicos do projeto PXD012998, focada na comparação entre o Carcinoma Seroso de Alto Grau (HGSC), o Carcinoma Endometrioide (EC) e o tecido da tuba uterina (Fallopian).

## 1. Fonte de Dados
Os dados originais foram obtidos do repositório PRIDE:
* **Projeto:** [PXD012998](https://www.ebi.ac.uk/pride/archive/projects/PXD012998)
* **Artigo de Referência:** Dieters-Castator et al. (2019).

### Amostras Utilizadas:
- **Fallopian Tissue (Controlo):** t1 a t4 (4 amostras)
- **HGSC:** t5 a t14 (10 amostras)
- **EC:** t15 a t24 (10 amostras)

## 2. Aquisição de Dados (Orbitrap Elite)
Embora o projeto original contenha dados de duas plataformas, esta reanálise foca-se exclusivamente nos dados obtidos pelo espectrómetro Orbitrap Elite.
```bash
wget [https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD012998/MaxQuant_proteinGroups_OrbitrapEliteData.txt](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD012998/MaxQuant_proteinGroups_OrbitrapEliteData.txt)

```

## 3. Processamento e Limpeza de Dados

A limpeza dos dados seguiu critérios rigorosos para garantir a qualidade das comparações estatísticas:

### Filtragem de Colunas:

* Manutenção das colunas de anotação: **Protein names**, **Protein ID** e **Gene names**.
* Manutenção exclusiva das intensidades **LFQ (Label-Free Quantification)**.
* **Tratamento de Duplicados:** Em amostras com múltiplas medições, deu-se prioridade à versão reprocessada (sufixo .2). As amostras de reinjeção (sufixo b) foram eliminadas.

### Filtragem de Linhas:

* **Contaminantes e Reversos:** Remoção de entradas que começam por "CON_" (contaminantes técnicos) e "REV__" (sequências reversas).
* **Múltiplas Identificações:** Eliminação de entradas com múltiplos Gene Names, exceto no caso de isoformas identificadas por hífen.

## 4. Estrutura de Ficheiros no Repositório

### Scripts (Linguagem R):

* `data_structure.R`: Organização e limpeza inicial das colunas do MaxQuant.
* `histogram_log2.R`: Script para visualização da distribuição das intensidades.
* `Treatment.R`: Processamento final das intensidades LFQ e gestão de réplicas.
* `DE_2.R`: Análise de expressão diferencial utilizando o pacote limma.
* `GO_Analysis.R`: Análise de enriquecimento funcional (Gene Ontology).

### Ficheiros e Resultados:

* `Sample_information.xlsx`: Metadados e anotação das amostras.
* `Treatment.zip`: Dados processados da plataforma Orbitrap Elite.
* `Normalization_2.zip`: Resultados da normalização e gráficos exploratórios (Boxplots e Volcano Plots).
* `GO_Analysis.zip`: Resultados das análises funcionais e dotplots de enriquecimento.

## 5. Metodologia Estatística

A principal distinção desta reanálise face ao estudo original reside na aplicação do **valor de p ajustado (Benjamini-Hochberg)**. Enquanto os autores originais basearam as suas conclusões em valores de p nominais, este projeto utiliza a correção para múltiplos testes para garantir a robustez dos biomarcadores e controlar a taxa de falsas descobertas (FDR).

---

**Autores:** Sara Cunha, Diogo Belbute, Luís Mota

**Data:** Dezembro 2025
