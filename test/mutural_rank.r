BiocManager::install(c("GenomicFeatures", "bigmemory"))
library(bigmemory) #necesary for first time mutual ranks calculation.
library(GenomicFeatures)

####normalization####
# "matrix" should be a count file like the output of htseq
matrix <- read.csv("samll_test.csv") 
txdbMp <- makeTxDbFromGFF("C:/Users/froma/Documents/Resources/Marchantia genome v5 2020/MpTak1v5.1_r1(new).gff") #this is the .gff file from the Marpolbase
df <- as.list(txdbMp)
df <- df$transcripts
Mp <- transcriptLengths(txdbMp)
genes <- rownames(matrix) #your gene list, usually is the same than unique(Mp$gene_id)
lengths <- Mp[match(genes,Mp$gene_id),]$tx_len #pick only the first isoform for each gene length
x <- matrix/lengths 
tpm <- t( t(x) * 1e6 / colSums(x) ) #counts to TPM
tpm <- log(tpm+0.125,2) #log-transformation
tpm <- tpm-rowMeans(tpm) #zero-centering

#### Calculate Pearson Correlation Coeficients ####
l <- dim(tpm)[1] ## length of the matrix
resultsm <- matrix(nrow=l,ncol=l) #generate an empty square matrix
n=1 ## for progress counting purposes only
m=1
#from here, the process will take a considerable amount of time.
for(i in 1:l){
  test <- apply(tpm,1,cor,tpm[i,])
  rank <- as.integer(rank(-test,na.last='keep')) ## replace corr values for ranks. Values are inverted so most coexpressed pair is 1 and less is l
  resultsm[i,] <- rank ## put ranks into data.frame
  #resultsm[i,] <- test ## modified to compute normal corr
  n=n+1
  if(n>m) {m=m+100 ## you can follow the progress using this counter
  message(paste('Progress ',round(n*100/l,2),'%...',sep=''))} 
}
#write.csv(resultsm,'~/Mp286+log.rank') ## not really necesary, but serves to save an intermediate output... takes a while

#resultsm<-read.csv('Mp274+log.rank',row.names = 1) #if the resultsm was pre-calculated
#resultsm<-as.matrix(resultsm)

#### Calculate Mutual Ranks ####
mutual <- big.matrix(nrow=l,ncol=l,type="double") # Big Marix is usually required to avoid a memory crush, but probably can be implemented in a normal matrix with a good computer
n=0 ## for progress counting purposes only
m=0
for(i in 1:l){
  mutual[i,] <- round((resultsm[i,]*resultsm[,i])^(1/2),1) ## here is the calculation as described in Obashashi et al. 2018 PCP.
  n=n+1
  if(n>m) {m=m+100 ## you can follow the progress using this
  message(paste('Progress ',round(n*100/l,2),'%....',sep=''))} 
}
write.big.matrix(mutual,'~/Mp_244+2.mutual') #this is the final matrix, the names of rows and columns are not stored and are the same as the original matrix.
#mutual <- read.big.matrix('~/R/cor_new/matrix_versions/Mp_244+.mutual',type='double') # if the mutual rank was pre-calculated

