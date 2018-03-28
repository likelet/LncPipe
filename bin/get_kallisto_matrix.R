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



#write reads count matrix of protein coding and lncRNA
#collapse expression in to gene level
row.names(map.file)<-map.file[,2]
map.file=map.file[pc.names,]

genelist<-unique(map.file$gene)
trans.kallisto.coun.out=apply(trans.kallisto.count,2,tapply,map.file$gene,sum)
row.names(trans.kallisto.coun.out)=levels(as.factor(genelist))

trans.kallisto.coun.out<-data.frame(ID=levels(as.factor(genelist)),data.frame(trans.kallisto.coun.out))
trans.kallisto.TPM.out=apply(trans.kallisto.TPM,2,tapply,map.file$gene,sum)
row.names(trans.kallisto.TPM.out)=levels(as.factor(genelist))
trans.kallisto.TPM.out<-data.frame(ID=levels(as.factor(genelist)),data.frame(trans.kallisto.TPM.out))

#annoted gene type
map.file2<-map.file[,-2]
names(map.file2)[1]="ID"
map.file2=unique(map.file2)
row.names(map.file2)=map.file2[,1]
typelist=map.file2[levels(as.factor(genelist)),2]
trans.kallisto.count.merge <- cbind(trans.kallisto.coun.out[,1],typelist,trans.kallisto.coun.out[,-1])
fwrite(trans.kallisto.count.merge ,file="kallisto.count.txt",sep="\t",quote=F,row.names=F)

#write tpm  matrix of  protein coding and lncRNA
trans.kallisto.TPM.merge <- cbind(trans.kallisto.coun.out[,1],typelist,trans.kallisto.TPM.out[,-1])
#annoted gene type
fwrite(trans.kallisto.TPM.merge ,file="kallisto.tpm.txt",sep="\t",quote=F,row.names=F)

