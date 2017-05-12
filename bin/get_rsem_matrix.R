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

pc.sample.order=c(grep("A",pc.samples),grep("D",pc.samples),grep("B",pc.samples),grep("C",pc.samples))
pc.rsem.fpkm=pc.rsem.fpkm[,pc.sample.order]
write.table(pc.rsem.fpkm,file="lncRNA/pc.rsem.fpkm.txt",sep="\t",quote=F)

pc.rsem.count=pc.rsem.count[,pc.sample.order]
write.table(pc.rsem.count,file="lncRNA/pc.rsem.count.txt",sep="\t",quote=F)

pc.rsem.TPM=pc.rsem.TPM[,pc.sample.order]
write.table(pc.rsem.TPM,file="lncRNA/pc.rsem.TPM.txt",sep="\t",quote=F)



lncRNA.sample.order=c(grep("A",lncRNA.samples),grep("D",lncRNA.samples),grep("B",lncRNA.samples),grep("C",lncRNA.samples))
lncRNA.rsem.fpkm=lncRNA.rsem.fpkm[,lncRNA.sample.order]
write.table(lncRNA.rsem.fpkm,file="lncRNA/lncRNA.rsem.fpkm.txt",sep="\t",quote=F)

lncRNA.rsem.count=lncRNA.rsem.count[,lncRNA.sample.order]
write.table(lncRNA.rsem.count,file="lncRNA/lncRNA.rsem.count.txt",sep="\t",quote=F)

lncRNA.rsem.TPM=lncRNA.rsem.TPM[,lncRNA.sample.order]
write.table(lncRNA.rsem.TPM,file="lncRNA/lncRNA.rsem.TPM.txt",sep="\t",quote=F)


