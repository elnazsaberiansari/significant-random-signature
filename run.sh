#!/bin/bash

## go to the directory where the code folder and data folder is there

cat $(pwd)/code/empiricalpvalue.R | R --slave --args $(pwd)/data/preprocessed-data.Rda $(pwd)/data/signatures.Rda $(pwd)/output

cat $(pwd)/code/randomSignatures.R | R --slave --args $(pwd)/data/preprocessed-data.Rda $(pwd)/data/signatures.Rda $(pwd)/output

cat $(pwd)/code/diffusion.R | R --slave --args $(pwd)/code/input/diffusion/genelist.txt $(pwd)/code/input/diffusion/genescore.txt $(pwd)/data/PPI.txt $(pwd)/output

cat $(pwd)/code/geneSubsetPvalue.R | R --slave --args $(pwd)/data/U133A_combat.h5 $(pwd)/code/input/geneSubsetPvalue/ACESGenesEntrezID.txt $(pwd)/code/input/geneSubsetPvalue/sigGenelist.txt $(pwd)/output
