# significant-random-signature


This repository contain the code for the paper :
> Significant random signature reveals new biomarker for breast cancer, by Saberi Ansari et al.

<h3>Descripotion:</h3>

```empiricalpvalue.R```
computes empirical pvalue for each signature and writes it in "empiricalpvalue" file in output directory.

Inputs:

- [NKI dataset (preprocessed-data.Rda)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240#s5)
- [dataset of 48 published signatures (signatures.Rda)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240#s5)
- directory of output


```randomSignatures.R```
compute random signatures for each of the 48 published signatures in NKI dataset and writes their nominal and empirical pvalues in "randomSignature_empiricalpvalue" file in the ouput folder. It scores the genes that  appear in significant random signatures and writes the genes and their scores in "significantRandomGenes" file in output folder.

Inputs:

- [NKI dataset (preprocessed-data.Rda)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240#s5)
- [dataset of 48 published signatures(signatures.Rda)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240#s5)
- directory of output


```diffusion.R```
Computes the diffusion score of the genes' score using [PPI network](https://we.tl/t-bLmPnAwgsa) and writes it in the "diffusionGeneScores.txt" in the output folder. Also computes the empirical p-values for diffusion scores and finds the significant genes which have an empirical p-value of less than 0.05, and writes them in the empiricalDS.txt and sigGenesDS.txt files in the output folder.

Inputs:

- ensemble id of the genes and their corresponding scores which is the output of randomSignatures.R
- PPI scores 
- directory of output


```geneSubsetPvalue.R```
Computes the nominal and empirical p-vlaues for the subset of genes and writes it in geneSetPvalue.txt file in output folder

Inputs:

- [ACES dataset](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3870302/)
- a file contatining EntrezID of ACESGenes
- EntrezID of the subset of genes which we want to check if they can seperate the poor and good prognosis of patients
- directory of output

