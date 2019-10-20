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
Computes the diffusion score of the genes' score using [PPI network](https://string-db.org) and writes it in the "diffusionGeneScores.txt" in the output folder. Also computes the empirical p-values for diffusion scores and finds the significant genes and writes them in the "empiricalDS.txt" and "sigGenesDS.txt" files in the output folder.

Inputs:

- files containing the ensemble id of the genes and their corresponding scores which are in represented in the output file of randomSignatures.R
- [PPI scores](https://string-db.org)
- directory of output


```geneSubsetPvalue.R```
Computes the nominal and empirical p-vlaues for the subset of genes and writes it in geneSetPvalue.txt file in output folder

Inputs:

- [ACES dataset (U133A_combat.h5)](https://ccb.nki.nl/software/aces/)
- a file contatining EntrezID of ACESGenes
- a file containing EntrezID of the subset of genes that we want to check if they can seperate the poor and good prognosis of patients
- directory of output



An additional dataset related to this analysis is also available at [IPM website](http://bs.ipm.ir/softwares/srs/). Ten breast cancer patients undergoing curative resection are included in this dataset. The expression of TAT (Tyrosine aminotransferase) in tumor, normal and cell lines are gathered in tat-2.xlsx file.
