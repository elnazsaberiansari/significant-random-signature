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
# output="/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019/code/output"

####################################################

getExp <- function(data, csig) {
  data$data <- data$data[(data$genes %in% csig$sig),, drop=F] 
  data$genes <- data$genes[(data$genes %in% csig$sig)]
  data$geneDescr <- data$geneDescr[(data$genes %in% csig$sig),, drop=F]
  data$idPlat <- NULL
  return(data)
}

####################################################
  
  pvaluecox <- numeric(length(cancer.signatures))
  empericalp <- numeric(length(cancer.signatures))
  for (signum in 1:length(cancer.signatures)) {
  	sigExp <- getExp(study.med, cancer.signatures[[signum]])
  	pc1 <- prcomp(t(sigExp$data))$x[,1]
  	med <- median(pc1)
    partition <- as.numeric(pc1 > med)
    pvaluecox[signum] <- log10(summary(coxph(sigExp$survival ~ partition))$logtest["pvalue"])
    surv1 <- sigExp$survival[partition == 1]
  	surv2 <- sigExp$survival[partition == 0]
  	
  	numofpermutation <- 1000
  	randompartition <- list()
  	pvaluecoxr <- numeric(numofpermutation)
  	
  	for(i in 1:numofpermutation) {
  		randompartition <- rep(FALSE, length(partition))
  		randomindex <- sample(1:length(partition), nrow(surv1))
  		randompartition[randomindex] <- TRUE
  		pvaluecoxr[i] <- log10(summary(coxph(sigExp$survival ~ randompartition))$logtest["pvalue"])	
  	}
  	empericalp[signum] = sum(abs(as.numeric(pvaluecoxr)) >= abs(pvaluecox[signum]))/numofpermutation
  }
  
  write.table(cbind(pvaluecox, empericalp), file = file.path(output,"empiricalpvalue"), row.names=TRUE, col.names=TRUE)


