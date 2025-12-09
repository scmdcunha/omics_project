# OMICS PROJECT

## Database em:

https://www.ebi.ac.uk/pride/archive/projects/PXD012998

## Samples_Information.xlsx

### Para obter o ficheiro:

https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD012998/Sample_information.xlsx

### Amostras:

- Fallopian Tissue: t1 - t4 (4 amostras)
- HGSC: t5 - t14 (10 amostras)
- EC: t15 - t24 (10 amostras)

## maxquant_OrbitrapElite_genes.txt e MaxQuant de QExactivePlusData.txt

### Download do MaxQuant_proteinGroups_OrbitrapEliteData.txt:
```bash
wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD012998/MaxQuant_proteinGroups_OrbitrapEliteData.txt
```

### Download do MaxQuant_proteinGroups_QExactivePlusData.txt:
```bash
wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD012998/MaxQuant_proteinGroups_QExactivePlusData.txt
```

### Colunas a manter:

- Protein names
- Protein ID
- Gene names
- LFQ intensity T**

### Linhas a eliminar:

https://github.com/scmdcunha/omics_project/issues/1#issue-3686795596

(Múltiplas identificações (em Gene Names) exceto isoformas (separadas por hífen); "CON_" E "REV__" (em Protein IDs))
