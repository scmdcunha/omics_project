# OMICS PROJECT

## Samples Information

### Para obter o ficheiro:

```bash
https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD012998/Sample_information.xlsx
```

### Amostras:

- Fallopian Tissue: t1 - t4
- HGSC: t5 - t14
- EC: t15 - t24

## maxquant_OrbitrapElite_genes.csv

### Para obter o ficheiro:

```bash
wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD012998/MaxQuant_proteinGroups_OrbitrapEliteData.txt
```

### Colunas a manter:

- Protein names
- Protein ID
- Gene names
- LFQ intensity T**

### Linhas a eliminar:

https://github.com/scmdcunha/omics_project/issues/1#issue-3686795596

(Múltiplas identificações (em Gene Names) exceto isoformas (separadas por hífen); "CON_" E "REV__" (em Protein IDs))
