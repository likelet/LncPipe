library(data.table)
kallisto.files <- list.files("./",pattern="*.tsv")
map.file <- read.table("./map.file",header=F,sep="\t")
names(map.file) <- c("ID","Type")

trans.kallisto.count <- c()
trans.kallisto.TPM <- c()
pc.names <- c()
for(pc in kallisto.files){
  pc.exp <- read.table(pc,header=T,sep="\t",row.names=1)
  trans.kallisto.count <- cbind(trans.kallisto.count,pc.exp[,3])
  trans.kallisto.TPM <- cbind(trans.kallisto.TPM,pc.exp[,4])
  pc.names <- rownames(pc.exp)
}
#assign gene names
rownames(trans.kallisto.count) <- pc.names
rownames(trans.kallisto.TPM) <- pc.names

#parsing samplenames
pc.samples <- unlist(lapply(kallisto.files,function(x){unlist(strsplit(x,"_"))[1]}))
colnames(trans.kallisto.count) <- pc.samples
colnames(trans.kallisto.TPM) <- pc.samples



#write reads count matrix of protein coding and lncRNA
trans.kallisto.count <- data.frame(ID=row.names(trans.kallisto.count),trans.kallisto.count)
#annoted transcript type
trans.kallisto.count.merge <- merge(map.file,trans.kallisto.count,by="ID")
fwrite(trans.kallisto.count.merge ,file="kallisto.count.txt",sep="\t",quote=F,row.names=F)
#write tpm  matrix of  protein coding and lncRNA
trans.kallisto.TPM <- data.frame(ID=row.names(trans.kallisto.TPM),trans.kallisto.TPM)
trans.kallisto.TPM.merge <- merge(map.file,trans.kallisto.TPM,by="ID")
#annoted transcript type
fwrite(trans.kallisto.TPM.merge ,file="kallisto.tpm.txt",sep="\t",quote=F,row.names=F)

