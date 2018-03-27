library(data.table)
kallisto.files <- list.files("./",pattern="*.tsv")
map.file <- read.table("./map.file",header=F,sep="\t")
names(map.file) <- c("gene","ID","Type")
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


trans.kallisto.count<-trans.kallisto.count[map.file$ID,]
trans.kallisto.TPM<-trans.kallisto.TPM[map.file$ID,]
#write reads count matrix of protein coding and lncRNA
#collapse expression in to gene level
trans.kallisto.coun.out=apply(trans.kallisto.count,2,tapply,map.file$gene,sum)
genelist<-unique(map.file$gene)
trans.kallisto.coun.out<-data.frame(ID=genelist,data.frame(trans.kallisto.coun.out))
trans.kallisto.TPM.out=apply(trans.kallisto.TPM,2,tapply,map.file$gene,sum)
trans.kallisto.TPM.out<-data.frame(ID=genelist,data.frame(trans.kallisto.TPM.out))

#annoted gene type
map.file2<-map.file[,-2]
names(map.file2)[1]="ID"
map.file2=map.file2[genelist,]
trans.kallisto.count.merge <- merge(map.file2,trans.kallisto.coun.out,by.x="ID")
fwrite(trans.kallisto.count.merge ,file="kallisto.count.txt",sep="\t",quote=F,row.names=F)

#write tpm  matrix of  protein coding and lncRNA
trans.kallisto.TPM.merge <- merge(map.file2,trans.kallisto.TPM.out,by="ID")
#annoted gene type
fwrite(trans.kallisto.TPM.merge ,file="kallisto.tpm.txt",sep="\t",quote=F,row.names=F)

