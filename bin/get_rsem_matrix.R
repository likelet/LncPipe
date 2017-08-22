library(data.table)
rsem.files=dir("RSEM")

pc.rsem=rsem.files[grep("pc.genes.results",rsem.files)]
lncRNA.rsem=rsem.files[grep("lncRNA.genes.results",rsem.files)]

pc.rsem.fpkm=c()
pc.rsem.count=c()
pc.rsem.TPM=c()
pc.names=c()
for(pc in pc.rsem){
  pc.exp=read.table(paste("RSEM/",pc,sep=""),header=T,sep="\t",row.names=1)
  pc.rsem.fpkm=cbind(pc.rsem.fpkm,pc.exp[,6])
  pc.rsem.count=cbind(pc.rsem.count,pc.exp[,4])
  pc.rsem.TPM=cbind(pc.rsem.TPM,pc.exp[,4])
  pc.names=rownames(pc.exp)
}

rownames(pc.rsem.fpkm)=pc.names
rownames(pc.rsem.count)=pc.names
rownames(pc.rsem.TPM)=pc.names

pc.samples=unlist(lapply(pc.rsem,function(x){unlist(strsplit(x,"_"))[1]}))
colnames(pc.rsem.fpkm)=pc.samples
colnames(pc.rsem.count)=pc.samples
colnames(pc.rsem.TPM)=pc.samples

lncRNA.rsem.fpkm=c()
lncRNA.rsem.count=c()
lncRNA.rsem.TPM=c()
lncRNA.names=c()
for(lncRNA in lncRNA.rsem){
  lncRNA.exp=read.table(paste("RSEM/",lncRNA,sep=""),header=T,sep="\t",row.names=1)
  lncRNA.rsem.fpkm=cbind(lncRNA.rsem.fpkm,lncRNA.exp[,6])
  lncRNA.rsem.count=cbind(lncRNA.rsem.count,lncRNA.exp[,4])
  lncRNA.rsem.TPM=cbind(lncRNA.rsem.TPM,lncRNA.exp[,4])
  lncRNA.names=rownames(lncRNA.exp)
}

rownames(lncRNA.rsem.fpkm)=lncRNA.names
rownames(lncRNA.rsem.count)=lncRNA.names
rownames(lncRNA.rsem.TPM)=lncRNA.names

lncRNA.samples=unlist(lapply(lncRNA.rsem,function(x){unlist(strsplit(x,"_"))[1]}))
colnames(lncRNA.rsem.fpkm)=lncRNA.samples
colnames(lncRNA.rsem.count)=lncRNA.samples
colnames(lncRNA.rsem.TPM)=lncRNA.samples

#write rpkm matrix of protein coding
pc.sample.order=c(grep("A",pc.samples),grep("D",pc.samples),grep("B",pc.samples),grep("C",pc.samples))
pc.rsem.fpkm=pc.rsem.fpkm[,pc.sample.order]
pc.rsem.fpkm=data.frame(ID=row.names(pc.rsem.fpkm),pc.rsem.fpkm)
fwrite(pc.rsem.fpkm,file="pc.rsem.fpkm.txt",sep="\t",quote=F,row.names=F)

#write reads count matrix of protein coding
pc.rsem.count=pc.rsem.count[,pc.sample.order]
pc.rsem.count=data.frame(ID=row.names(pc.rsem.count),pc.rsem.count)
fwrite(pc.rsem.count,file="pc.rsem.count.txt",sep="\t",quote=F,row.names=F)
#write tpm matrix of protein coding
pc.rsem.TPM=pc.rsem.TPM[,pc.sample.order]
pc.rsem.TPM=data.frame(ID=row.names(pc.rsem.TPM),pc.rsem.TPM)
fwrite(pc.rsem.TPM,file="pc.rsem.TPM.txt",sep="\t",quote=F,row.names=F)


#write rpkm matrix of long non coding RNAS
lncRNA.sample.order=c(grep("A",lncRNA.samples),grep("D",lncRNA.samples),grep("B",lncRNA.samples),grep("C",lncRNA.samples))
lncRNA.rsem.fpkm=lncRNA.rsem.fpkm[,lncRNA.sample.order]
lncRNA.rsem.fpkm=data.frame(ID=row.names(lncRNA.rsem.fpkm),lncRNA.rsem.fpkm)
fwrite(lncRNA.rsem.fpkm,file="lncRNA.rsem.fpkm.txt",sep="\t",quote=F,row.names=F)
#write reads count matrix of  long non coding RNAS
lncRNA.rsem.count=lncRNA.rsem.count[,lncRNA.sample.order]
lncRNA.rsem.count=data.frame(ID=row.names(lncRNA.rsem.count),lncRNA.rsem.count)
fwrite(lncRNA.rsem.count,file="lncRNA.rsem.count.txt",sep="\t",quote=F,row.names=F)
#write tpm  matrix of long non-coding RNAS
lncRNA.rsem.TPM=lncRNA.rsem.TPM[,lncRNA.sample.order]
lncRNA.rsem.TPM=data.frame(ID=row.names(lncRNA.rsem.TPM),lncRNA.rsem.TPM)
fwrite(lncRNA.rsem.TPM,file="lncRNA.rsem.TPM.txt",sep="\t",quote=F,row.names=F)


