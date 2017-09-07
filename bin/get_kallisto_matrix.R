library(data.table)
kallisto.files <- list.files("./")

pc.kallisto <- kallisto.files[grep("kallisto.tsv",kallisto.files)]
lncRNA.kallisto <- kallisto.files[grep("kallisto.tsv",kallisto.files)]

pc.kallisto.count <- c()
pc.kallisto.TPM <- c()
pc.names <- c()
for(pc in pc.kallisto){
  pc.exp <- read.table(pc,header=T,sep="\t",row.names=1)
  pc.kallisto.count <- cbind(pc.kallisto.count,pc.exp[,3])
  pc.kallisto.TPM <- cbind(pc.kallisto.TPM,pc.exp[,4])
  pc.names <- rownames(pc.exp)
}
#assign gene names
rownames(pc.kallisto.count) <- pc.names
rownames(pc.kallisto.TPM) <- pc.names

#parsing samplenames
pc.samples <- unlist(lapply(pc.kallisto,function(x){unlist(strsplit(x,"_"))[1]}))
colnames(pc.kallisto.count) <- pc.samples
colnames(pc.kallisto.TPM) <- pc.samples

lncRNA.kallisto.count <- c()
lncRNA.kallisto.TPM <- c()
lncRNA.names <- c()
for(lncRNA in lncRNA.kallisto){
  lncRNA.exp <- read.table(lncRNA,header=T,sep="\t",row.names=1)
  lncRNA.kallisto.count <- cbind(lncRNA.kallisto.count,lncRNA.exp[,3])
  lncRNA.kallisto.TPM <- cbind(lncRNA.kallisto.TPM,lncRNA.exp[,4])
  lncRNA.names <- rownames(lncRNA.exp)
}

rownames(lncRNA.kallisto.count) <- lncRNA.names
rownames(lncRNA.kallisto.TPM) <- lncRNA.names

lncRNA.samples <- unlist(lapply(lncRNA.kallisto,function(x){unlist(strsplit(x,"_"))[1]}))
#add samplenames
colnames(lncRNA.kallisto.count) <- lncRNA.samples
colnames(lncRNA.kallisto.TPM) <- lncRNA.samples

#write reads count matrix of protein coding and lncRNA
pc.kallisto.count <- data.frame(ID=row.names(pc.kallisto.count),Type=rep("Coding",nrow(pc.kallisto.count)),pc.kallisto.count)
lncRNA.kallisto.count <- data.frame(ID=row.names(lncRNA.kallisto.count),Type=rep("lncRNA",nrow(pc.kallisto.count)),lncRNA.kallisto.count)
count.matrix <- rbind(pc.kallisto.count,lncRNA.kallisto.count)
fwrite(matrix,file="kallisto.count.txt",sep="\t",quote=F,row.names=F)
#write tpm  matrix of  protein coding and lncRNA
pc.kallisto.TPM <- data.frame(ID=row.names(pc.kallisto.TPM),Type=rep("lncRNA",nrow(pc.kallisto.count)),pc.kallisto.TPM)
lncRNA.kallisto.TPM <- data.frame(ID=row.names(lncRNA.kallisto.TPM),Type=rep("lncRNA",nrow(pc.kallisto.count)),lncRNA.kallisto.TPM)
TPM.matrix <- rbind(pc.kallisto.TPM,lncRNA.kallisto.TPM)
fwrite(TPM.matrix ,file="kallisto.tpm.txt",sep="\t",quote=F,row.names=F)

