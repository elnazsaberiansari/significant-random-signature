library(rhdf5)
library(survival)

## call from consule
args <- commandArgs()
path=args[4]
ds=args[5]
geneID=args[6]
genesubset=args[7]
output=args[8]

## call from R
# path="/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019"
# ds="data/U133A_combat.h5"
# geneID="/code/input/geneSubsetPvalue/ACESGenesEntrezID.txt"
# genesubset="/code/input/geneSubsetPvalue/sigGenelist.txt"
# output="/run/media/esaberia/Elnaz/allprojects/SignificantRandomSignature/SignificantRandomSignature_2019/code/output"

###################################################
####################################################

	# h5ls(file.path(path,"data/U133A_combat.h5"))
	numofpermutation <- 1000
	
  for(d in 1:2){
		# allDatasets = h5read("data/U133A_combat.h5", "/")
		if(d==1){
			Dataset = h5read(file.path(path,ds), "/U133A_combat_DMFS")
		}
		if(d==2){
			Dataset = h5read(file.path(path,ds), "/U133A_combat_RFS")
		}
	
		goodPatientDataset <- which(Dataset$PatientClassLabels==FALSE)
		goodPatientNum <- length(goodPatientDataset)
		poorPatientDataset <- which(Dataset$PatientClassLabels==TRUE)
		poorPatientNum <- length(poorPatientDataset)
		GoodPatientsExp <- Dataset$ExpressionData[,Dataset$PatientClassLabels==FALSE]
		PoorPatientsExp <- Dataset$ExpressionData[,Dataset$PatientClassLabels==TRUE]
		Genenum <- nrow(Dataset$GeneLabels)
	
		ACESENTREZID <- read.table(file.path(path,geneID) , header=TRUE)
		ACESENTREZID <- ACESENTREZID$ACESEntrezID
	
		SP <- readLines(paste(file.path(path,genesubset)))
		SP <- strsplit(SP,"\t")
		SP <- lapply(SP,as.numeric)
			
			for(pnum in 1:length(SP)) {
				pathwayName <- pnum
				SPGE <- SP[[pnum]]
			
				#ind is index of genes of the subpathway in the ACES database
				ind <- match(SPGE, ACESENTREZID)
				ind <- ind[!(ind %in% NA)]
			
				if(length(ind) > 1){
					subpathGoodPatientsExp <- GoodPatientsExp[ind,]
					subpathPoorPatientsExp <- PoorPatientsExp[ind,]
				
					#mean of each gene in Poor or Good
					meanSGP <- sapply(1:length(ind),function(x){mean(as.numeric(subpathGoodPatientsExp[x,]))})
					meanSPP <- sapply(1:length(ind),function(x){mean(as.numeric(subpathPoorPatientsExp[x,]))})
					
					subpathGoodPatientsExpm2 <- sapply(1:length(ind),function(x){abs(subpathPoorPatientsExp[x,]-meanSGP[x])^2})
					subpathPoorPatientsExpm2 <- sapply(1:length(ind),function(x){abs(subpathPoorPatientsExp[x,]-meanSPP[x])^2})
					
					tsubpathGoodPatientsExpm2 <- t(subpathGoodPatientsExpm2)
					tsubpathPoorPatientsExpm2 <- t(subpathPoorPatientsExpm2)
					
					SPGO2 <- sapply(1:length(ind),function(x){abs(subpathPoorPatientsExp[x,]* tsubpathGoodPatientsExpm2[x,])})
					SPPO2 <- sapply(1:length(ind),function(x){abs(subpathPoorPatientsExp[x,]* tsubpathPoorPatientsExpm2[x,])})
					
					Score1O2 <- sapply(1:nrow(SPGO2),function(x){sum(as.numeric(SPGO2[x,]))})
					Score2O2 <- sapply(1:nrow(SPPO2),function(x){sum(as.numeric(SPPO2[x,]))})
					
					nominalp <- t.test(Score1O2, Score2O2, var.equal=FALSE)$p.value
					nominalt <- t.test(Score1O2, Score2O2, var.equal=FALSE)$statistic

					totaldata <- cbind(subpathGoodPatientsExp, subpathPoorPatientsExp)
				
					empericaltr <- list()
					for(i in 1:numofpermutation) {
					  randomdata <- rep(FALSE, ncol(totaldata))
					  randomindex <- sample(1:ncol(totaldata), ncol(subpathGoodPatientsExp))
					  randomdata[randomindex] <- TRUE
					  subpathGoodPatientsExpr <- totaldata[ ,randomdata]
					  subpathPoorPatientsExpr <- totaldata[ ,!randomdata]
					  
					  meanSGPr <- sapply(1:length(ind),function(x){mean(as.numeric(subpathGoodPatientsExpr[x,]))})
					  meanSPPr <- sapply(1:length(ind),function(x){mean(as.numeric(subpathPoorPatientsExpr[x,]))})
					  
					  subpathGoodPatientsExpmr2 <- sapply(1:length(ind),function(x){ abs(subpathPoorPatientsExpr[x,]-meanSGPr[x])^2})
					  subpathPoorPatientsExpmr2 <- sapply(1:length(ind),function(x){ abs(subpathPoorPatientsExpr[x,]-meanSPPr[x])^2})
					  
					  tsubpathGoodPatientsExpmr2 <- t(subpathGoodPatientsExpmr2)
					  tsubpathPoorPatientsExpmr2 <- t(subpathPoorPatientsExpmr2)
					  
					  SPGOr2 <- sapply(1:length(ind),function(x){abs(subpathPoorPatientsExpr[x,]* tsubpathGoodPatientsExpmr2[x,])})
					  SPPOr2 <- sapply(1:length(ind),function(x){abs(subpathPoorPatientsExpr[x,]* tsubpathPoorPatientsExpmr2[x,])})
					  
					  Score1Or2 <- sapply(1:nrow(SPGOr2),function(x){sum(as.numeric(SPGOr2[x,]))})
					  Score2Or2 <- sapply(1:nrow(SPPOr2),function(x){sum(as.numeric(SPPOr2[x,]))})
					  
					  empericaltr <- cbind(empericaltr, t.test(Score1Or2, Score2Or2, var.equal=FALSE)$statistic)
					}
					
					empericalp = sum(abs(as.numeric(empericaltr)) >= abs(as.numeric(nominalt))) / numofpermutation
					write.table(cbind(nominalp, empericalp), file=file.path(output,paste0("geneSetPvalue.txt")), append=TRUE, row.names=FALSE, col.names= FALSE)
					
				}
		}
	}
