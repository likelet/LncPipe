basic_charac=read.table("lncRNA/basic_charac.txt",header=F,sep="\t")

boxplot(basic_charac[which(basic_charac[,3]=="known"),4],basic_charac[which(basic_charac[,3]=="novel"),4],basic_charac[which(basic_charac[,3]=="protein_coding"),4])

known.cp=ecdf(c(basic_charac[which(basic_charac[,3]=="known"),4],1))
novel.cp=ecdf(c(basic_charac[which(basic_charac[,3]=="novel"),4],1))
pc.cp=ecdf(basic_charac[which(basic_charac[,3]=="protein_coding"),4])
pdf("lncRNA/compare_coding_potential.pdf")
plot(pc.cp,lwd=4,col="red",xlim=c(0,1),main="Coding Potential",ylab="CDF",xlab="Coding Probablity (CPAT)")
lines(known.cp,lwd=4,col="orange",xlim=c(0,1))
lines(novel.cp,lwd=4,col="blue",xlim=c(0,1))
dev.off()

novel.size=density(basic_charac[which(basic_charac[,3]=="novel"),5])
known.size=density(basic_charac[which(basic_charac[,3]=="known"),5])
pc.size=density(basic_charac[which(basic_charac[,3]=="protein_coding"),5])

pdf("lncRNA/compare_transcript_length.pdf")
plot(known.size,col="orange",lwd=4,main="Transcript Length Distribution",xlab="Length",ylab="Density")
lines(novel.size,col="blue",lwd=4)
lines(pc.size,col="red",lwd=4)
dev.off()

novel.ec=table(basic_charac[which(basic_charac[,3]=="novel"),6])
novel.ec=novel.ec/sum(novel.ec)
known.ec=table(basic_charac[which(basic_charac[,3]=="known"),6])
known.ec=known.ec/sum(known.ec)
pc.ec=table(basic_charac[which(basic_charac[,3]=="protein_coding"),6])
pc.ec=pc.ec/sum(pc.ec)

ec=rbind(novel.ec[as.character(2:20)],known.ec[as.character(2:20)],pc.ec[as.character(2:20)])
colnames(ec)=2:20

pdf("lncRNA/compare_exon_number.pdf")
barplot(ec,beside=T,col=c("blue","orange","red"),xlab="Exon Count",ylab="Percentage of Transcripts")
dev.off()

