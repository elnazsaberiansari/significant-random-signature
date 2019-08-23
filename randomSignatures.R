library(survival)

## call from consule
args <- commandArgs()
load(args[4])
load(args[5])
output=args[6]

## call from R
# path="/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019"
# load(file.path(path,"/data/preprocessed-data.Rda"))
# load(file.path(path,"/data/signatures.Rda"))
# output="/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019/output"

###################################################
####################################################

GetRandGenes <- function(data, n) {
  ind <- sample(1:nrow(data$data), n)
  data$data <- data$data[ind,, drop=F]
  data$genes <- data$genes[ind]
  data$geneDescr <- data$geneDescr[ind,, drop=F]
  data$idPlat <- NULL
  return(data)
}

####################################################

getExp <- function(data, csig) {
  data$data <- data$data[(data$genes %in% csig$sig),, drop=F] 
  data$genes <- data$genes[(data$genes %in% csig$sig)]
  data$geneDescr <- data$geneDescr[(data$genes %in% csig$sig),, drop=F]
  data$idPlat <- NULL
  return(data)
}

####################################################

  significanrandGenes <- rep(0, length(study.med$genes))
  thr <- -10
  alpha <- 0
  pernum <- 1000
	for(signum in 1:length(cancer.signatures)){
	  for(r in 1:pernum){
    	  randsigExp <- GetRandGenes(study.med, length(cancer.signatures[[signum]]$sig))
    	  randpc1 <- prcomp(t(randsigExp $data))$x[,1]
    	  randmed <- median(randpc1)
    	  randpartition <- as.numeric(randpc1 > randmed)
    	  randpvaluecox <- log10(summary(coxph(randsigExp $survival ~ randpartition))$logtest["pvalue"])
    	  randsurv1 <- randsigExp$survival[randpartition == 1]
    	  randsurv2 <- randsigExp$survival[randpartition == 0]
    	  
    	  numofpermutation <- 1000
    	  randompartition <- list()
    	  pvaluecoxr <- numeric(numofpermutation)
    	  
    	  for(i in 1:numofpermutation) {
    	    rpartition <- rep(FALSE, length(randpartition))
    	    rindex <- sample(1:length(randpartition), nrow(randsurv1))
    	    rpartition[rindex] <- TRUE
    	    pvaluecoxr[i] <- log10(summary(coxph(randsigExp$survival ~ rpartition))$logtest["pvalue"])	
    	  }
    	  empericalp = sum(abs(as.numeric(pvaluecoxr)) >= abs(randpvaluecox))/numofpermutation
    	  write.table(cbind(paste0("r",signum,"_",r),randpvaluecox, empericalp), file = file.path(output,"randomSignature_empiricalpvalue"), row.names=FALSE, col.names=FALSE, append = TRUE)
    	  
    	    ## computing number of genes seen in significqnt rqndom signqture
      	  if (pvaluecoxr < thr && pvaluecoxr == alpha) {
      	    genesindex <- study.med$genes %in% randsigExp$genes
      	    significanrandGenes[genesindex] <- significanrandGenes[genesindex] + 1 
      	  }
	  }
	}
  write.table(cbind(study.med$geneDescr$geneId, study.med$geneDescr$symb, significanrandGenes), file = file.path(output,"significantRandomGenes"), row.names=FALSE, col.names=FALSE)
  