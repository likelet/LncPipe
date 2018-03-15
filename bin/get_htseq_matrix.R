library(data.table)
htseq.files <- list.files("./",pattern="*.htseq.count")
map.file <- read.table("./map.file",header=F,sep="\t")
names(map.file) <- c("ID","Type")

trans.htseq.count <- c()

pc.names <- c()
for(pc in htseq.files){
  pc.exp <- read.table(pc,header=T,sep="\t",row.names=1)
  trans.htseq.count <- cbind(trans.htseq.count,pc.exp[,1])
  pc.names <- rownames(pc.exp)
}
#assign gene names
rownames(trans.htseq.count) <- pc.names

#parsing samplenames
pc.samples <- unlist(lapply(htseq.files,function(x){unlist(strsplit(x,"_"))[1]}))
colnames(trans.htseq.count) <- pc.samples
#remove the useless line
rans.htseq.count<-rans.htseq.count[1:(nrow(rans.htseq.count)-5),]


#write reads count matrix of protein coding and lncRNA
trans.htseq.count <- data.frame(ID=row.names(trans.htseq.count),trans.htseq.count)
#annoted transcript type
trans.htseq.count.merge <- merge(map.file,trans.htseq.count,by="ID")
fwrite(trans.htseq.count.merge ,file="htseq.count.txt",sep="\t",quote=F,row.names=F)


