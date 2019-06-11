### Generate toy data

getwd()
# under the SNP_VarianceHet_Tools_withReedik directory



# get pheno for toy data
samplefile1 <- read.table("./Stage1TestingFiles/Study1/BMI_1.txt", head=T)
# provide site specific sample size
n_vect <- c(1000, 500, 3000, 1000, 1000)

# create pheno data for each new study 2-6
for (i in 1:5){
	
	samplefile <- data.frame("FID" = rep(0, n_vect[i]), 
							 "IID" = 1:n_vect[i], 
							 "bmi_lg" = sample(samplefile1$bmi_lg, n_vect[i], replace=T))
							 
	write.table(samplefile, file = paste("./Stage1TestingFiles/Study", 1+i, "/BMI_", 1+i, ".txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
	
}



# get tfam for toy data
sampleTpedfam1 <- read.table("./Stage1TestingFiles/Study1/tped_test1.tfam", head=F)

# create pheno data for each new study 2-6
for (i in 1:5){
	
	sampleTpedfam <- data.frame("FID" = rep(0, n_vect[i]), 
								"IID" = 1:n_vect[i], 
								"M" = rep(0, n_vect[i]), 
								"F" = rep(0, n_vect[i]), 
								"Gender" = rbinom(n_vect[i],2,0.4), 
								"Pheno"= rep(-9, n_vect[i]))
								
	write.table(sampleTpedfam, file = paste("./Stage1TestingFiles/Study", 1+i, "/tped_test", 1+i, ".tfam", sep=""), row.names=F,  col.names=F, sep="\t")
	
}




# get tfam for toy data
sampleTped1 <- read.table("./Stage1TestingFiles/Study1/marker_1.tped", head=F)

# create pheno data for each new study 2-6
for (i in 1:5){

	if (n_vect[i] < dim(sampleTped1)[2]-4){
	
		sampleTped <- sampleTped1[,c(1:4, sample(5:dim(sampleTped1)[2], n_vect[i], replace=T))]	
			
	} else {

		sampleTped[,1:4] <- sampleTped1[,1:4]
		sampleTped[,5:(n_vect[i]+4)] <- rbinom(n_vect[i], 2,  0.3)									
	}
	
	write.table(sampleTped, file = paste("./Stage1TestingFiles/Study", 1+i, "/marker_", 1+i, ".tped", sep=""), quote=F, row.names=F,  col.names=F)
	
}

