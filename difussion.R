
Diffusion <- function(GGW, Beta) {
	######## Construct Diffusion Kernel ########
  	D <- diag(colSums(GGW))
  	L<-GGW-D
  	EG <- eigen(L, sym="TRUE")
  	E <- exp(Beta * EG$values)
  	K <- (EG$vectors %*% diag(E) %*% solve(EG$vectors))
  	K[which(K<0)] <- 0
    return(K)
}

####################################################

ComputeSigGenesByDiffusion <- function (GGW, GeneList, GeneScores, Beta, Alpha) {

  	######## Computed Diffused Score with function Diffusion ########
  	K <- Diffusion(GGW, Beta)
  	DS <- GS$V1 %*% K
  	
  	write.table(K, file=file.path(output,paste0("diffusionKernel.txt")), row.names=FALSE, col.names=FALSE)
  	write.table(cbind(GS, t(DS)), file = file.path(output,paste0("diffusionGeneScores.txt")), row.names=FALSE, col.names=TRUE)
  	# K <- read.table(sep="\t", file=file.path(output,paste0("diffusionKernel-beta=", Beta,".txt")), header = FALSE)
  	# score <- read.table(sep="\t", file=file.path(output,"GGWADj"), header = TRUE)
  	
  
  	######## Identify Significant Genes ########
  	nPerm <- 1000
  	rDS <- numeric(geneNum)
  	empericalDS <- list()
  	for(i in 1: nPerm) {
  		rDS <- sample(GS) %*% K
  		empericalDS <- cbind(empericalDS, t(rDS))
  	}
	
	empericalp = numeric(geneNum)
	for(i in 1: geneNum){
		empiricalp[i] = sum(abs(as.numeric(empericalDS[i,])) >= abs(as.numeric(DS[i])))/nPerm
	}
	
	sigGenes <- which( empiricalp <= Alpha)
	
	write.table(empiricalp, file=file.path(output,paste0("empiricalDS.txt")), row.names=FALSE, col.names=FALSE)
	write.table(sigGenes, file = file.path(output,paste0("sigGenesDS.txt")), row.names=FALSE, col.names=TRUE)
}

####################################################
  ## call from consule
  args <- commandArgs()
  GeneList=args[4]
  GeneScore=args[5]
  GeneGeneWeights=args[6]
  output=args[7]

  # GeneList <- "/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019/code/input/diffusion/genelist.txt"
  # GeneScore <- "/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019/code/input/diffusion/genescore.txt"
  # GeneGeneWeights <- "/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019/data/PPI.txt"
  # output <- "/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/sigtest/SignificantRandomSignature_files/output"
  
  GL <- read.table(GeneList, header=TRUE)
  GS <- read.table(GeneScore)
  PPI <- read.table(GeneGeneWeights, header=TRUE)
  
  ## check if all the GL are in PPI proteins
  # ind <- match(GL$genesID, PPI$protein2)
  # indna <- ind %in% NA
  # GList <- GL[!indna,]
  # GScore <- GS[!indna,]
	 
	geneNum <- nrow(GL)
	GGWADj <- matrix(0, geneNum,geneNum)
	 
	 for (i in 1:geneNum) {
	  print(i)
		ind1 <- which(PPI$protein1 %in% GL$genesID[i]) # ind1 is the indices of GL$genesID[i] in PPI$protein1
	 	ind2 <- match(PPI$protein2[ind1],GL$genesID)
	 	# ind2 <- which(PPI$protein2[ind1] %in% GL$genesID)
	 	indna <- ind2 %in% NA
	 	ind3 <- ind2[!indna]
	 	ind4 <- ind1[!indna]
	 	GGWADj[i,ind3] <- (PPI$score[ind4])/1000
	 }
	 
   write.table(GGWADj, file=file.path(output,"GGWADj"), row.names=FALSE, col.names=FALSE)
   # GGWADj <- read.table(sep="\t", file=file.path(output,"GGWADj"), header = FALSE)
   
  	 
   Beta <- 3
   Beta <- Beta/100
   Alpha <- 0.05
   ComputeSigGenesByDiffusion(GGW=GGWADj, GeneList=GL, GeneScores=GS, Beta, Alpha)
