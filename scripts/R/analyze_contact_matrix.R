# R --slave --args file bin filedistances filecolumns projectdir < analyze_contact_matrix.R

library(Heatplus)

library(gplots)

library(RColorBrewer)

library(MCMCpack)

orgPal<-brewer.pal(9,"YlGnBu")

args = commandArgs(TRUE)





filename 	<- 	args[1]



bin 		<- 	args[2]



filedistances 	<- 	args[3]



filecolumns 	<-	args[4]



filenamecvg 	<-	args[5]





projectdir	<-	args[6]





setwd(projectdir)





t<-as.matrix(read.table(paste("matrices",filecolumns,sep="/")))





# read contact matrix

data_matrix<-data.matrix(read.table(paste("matrices",filename,sep="/"),header=TRUE,row.names= 1))

data_matrix_dim<-dim(data_matrix)[1]





# read cvg matrix

filenamecvg<-paste(filename,"cvg",sep = '.')

data_matrixcvg<-data.matrix(read.table(paste("matrices",filenamecvg,sep="/"),header=TRUE,row.names= 1))

data_matrixcvg_dim<-dim(data_matrixcvg)[1]



# read distances

distances<-abs(as.matrix(read.table(paste("matrices",filedistances,sep="/"))))



# read statistics

setwd("statistics")

read.stat<-list.files()

read1<-read.table(read.stat[1],comment.char = "", fill=TRUE)

read2<-read.table(read.stat[2],comment.char = "", fill=TRUE)

# generate bar-plots on reads statistics



#    tr1<-as.numeric(as.matrix(read1[1,4]))/1000000

#    ra1<-as.numeric(as.matrix(read1[4,2]))/1000000

#    na1<-as.numeric(as.matrix(read1[3,7]))/1000000



#    tr2<-as.numeric(as.matrix(read2[1,4]))/1000000

#    ra2<-as.numeric(as.matrix(read2[4,2]))/1000000

#    na2<-as.numeric(as.matrix(read2[3,7]))/1000000



#    intrachrom<-sum(diag(data_matrix))/1000000

#    interchrom<-sum(sum(data_matrix)-sum(diag(data_matrix)))/(2*1000000)



#    png('read.statistics.png',height = 500)

#    barplot(c(tr1+tr2,na1+na2,tr1,tr2,intrachrom+interchrom,intrachrom,interchrom),names.arg=c("Total","Not Aligning","Aligning r1","Aligning r2","Alignable pairs","intrachrom","interchrom"),ylab ="Number of reads (milions)", las=2, cex.names=0.8)

#    dev.off()







# generate plots on reads distances

breaks<-dim(data_matrix)[1]

png(paste(filedistances,'.hist.observed.png',sep=''))

hist(log(as.matrix(distances)), breaks = breaks,main = "Histogram of reads distances", xlab = "Log(distances bp)")

dev.off()



setwd("../")



# create directory for current bin 

bindir<-paste("bin",bin,sep=".")

dir.create(bindir)

setwd(bindir)



# create heatmap of observed interactions (split)

dir.create("observed")

setwd("observed")





p<-1



while (p<=dim(t)[1])



{



q<-1



while



(q<=dim(t)[1])



{



png(paste(t[p,1],"-",t[q,1],'.',bin,'.heatmap.observed.png',sep=''), width = 800, height = 800, units = "px", pointsize = 12)



heatmap.2 (log2(1+data_matrix[t[p,2]:t[p,3],t[q,2]:t[q,3]]), Colv=NA, Rowv=NA, dendrogram=c("none"),trace=c("none"), col=orgPal, ylab = t[p,1], xlab =t[q,1], labRow = NA, labCol =NA)



dev.off()



write.table(data_matrix[t[p,2]:t[p,3],t[q,2]:t[q,3]],file = paste(t[p,1],"-",t[q,1],'.',bin,'.heatmap.observed.txt',sep=''),sep = "\t")



q<-q+1



}



p<-p+1



}



setwd('../')



# Generation of Expected Matrix

orgPal<-brewer.pal(9,"BuPu")

dir.create("expected")

setwd("expected")



p<-1



while (p<=dim(t)[1])



{



data_matrix_chr<-(data_matrix[t[p,2]:t[p,3],t[p,2]:t[p,3]])



data_matrix_chr_dim<-dim(data_matrix_chr)[1]

x<-1



y=4.7*(x^(-0.7))



while (x < data_matrix_chr_dim) {



     ya=4.7*(x^(-0.7))



     y<-c(y,ya)



     x<-x+1

}



z<-y/c(data_matrix_chr_dim:1)



x<-length(y)



while (x > 1) {



     z<-c(z,z[1:(x-1)])



     x<-x-1



}



expected<-xpnd(z,data_matrix_chr_dim)



expected.f<-expected*10000000/sum(expected)



rownames(expected.f) <- rownames(data_matrix_chr)

colnames(expected.f) <- colnames(data_matrix_chr)

write.table(expected.f,file = paste(t[p,1],"-",t[p,1],'.',bin,'.expected.txt',sep=''),sep = "\t")



# create heatmap of expected-interactions (split)



png(paste(t[p,1],"-",t[p,1],'.',bin,'.heatmap.expected.png',sep=''), width = 800, height = 800, units = "px", pointsize = 12)



heatmap.2 (log2(1+(expected.f)), Colv=NA,Rowv=NA, dendrogram=c("none"),trace=c("none"), col=orgPal, ylab = t[p,1], xlab =t[p,1], labRow = NA, labCol = NA)



dev.off()



p<-p+1



}

setwd('../')





# Generation of Observed vs Expected Matrix

orgPal<-brewer.pal(9,"YlGn")

dir.create("observed.vs.expected")

setwd("observed.vs.expected")



p<-1



while (p<=dim(t)[1])



{



data_matrix_chr<-(data_matrix[t[p,2]:t[p,3],t[p,2]:t[p,3]])



data_matrix_dim<-dim(data_matrix_chr)[1]

x<-1



y=4.7*(x^(-0.7))



while (x < data_matrix_dim) {



     ya=4.7*(x^(-0.7))



     y<-c(y,ya)



     x<-x+1

}



z<-y/c(data_matrix_dim:1)



x<-length(y)



while (x > 1) {



     z<-c(z,z[1:(x-1)])



     x<-x-1



}



expected<-xpnd(z,data_matrix_dim)



data_matrix_chr.f<-data_matrix_chr*10000000/sum(data_matrix_chr)



expected.f<-expected*10000000/sum(expected)

obs.vs.exp<-(data_matrix_chr.f/expected.f)



rownames(obs.vs.exp) <- rownames(data_matrix_chr)

colnames(obs.vs.exp) <- colnames(data_matrix_chr)

write.table(obs.vs.exp,file = paste(t[p,1],"-",t[p,1],'.',bin,'.obs.vs.exp.txt',sep=''),sep = "\t")

# create heatmap of intra-interactions (split)

png(paste(t[p,1],"-",t[p,1],'.',bin,'.heatmap.intra-interaction.png',sep=''), width = 800, height = 800, units = "px", pointsize = 12)



heatmap.2 (log2(1+(obs.vs.exp)), Colv=NA,Rowv=NA, dendrogram=c("none"),trace=c("none"), col=orgPal, ylab = t[p,1], xlab =t[p,1], labRow = NA, labCol = NA)



dev.off()



p<-p+1



}

setwd('../')







# Correction for Coverage



data_matrixcvg<-data_matrixcvg/sum(data_matrix)



orgPal<-brewer.pal(9,"Reds")

dir.create("coverage-correction")

setwd("coverage-correction")



p<-1



while (p<=dim(t)[1])



{



q<-1



while



(q<=dim(t)[1])



{



png(paste(t[p,1],"-",t[q,1],'.',bin,'.heatmap.cvg.png',sep=''), width = 800, height = 800, units = "px", pointsize = 12)



cvg<-data_matrixcvg[t[p,2]:t[p,3],t[q,2]:t[q,3]]



heatmap.2 (log2(1+(cvg)), Colv=NA, Rowv=NA, dendrogram=c("none"),trace=c("none"), col=orgPal, ylab = t[p,1], xlab =t[q,1], labRow = NA, labCol =NA)



dev.off()



write.table(cvg,file = paste(t[p,1],"-",t[q,1],'.',bin,'.heatmap.cvg.txt',sep=''),sep = "\t")



q<-q+1



}



p<-p+1



}



setwd('../')



# Generation of Correlation Matrix

orgPal<-brewer.pal(9,"OrRd")

dir.create("correlation")

setwd("correlation")

x<- sqrt(data_matrix)

x[is.na(x)] <- 0

i<-1

matcor <- matrix(data = NA, nrow = dim(x)[1], ncol = dim(x)[2])

while (i <= dim(x)[1]) {

	j<-1

while (j <= dim(x)[1]) 

{

	 matcor[i,j] <- try(cor(x[i,],x[,j], use = "complete.obs"), silent = TRUE)

	j<-j+1

}

	i<-i+1
}



matcor<-matrix(matcor[2:(length(matcor))],dim(x)[1],dim(x)[1])







# create heatmap of correlations (split)



p<-1

while (p<=dim(t)[1])


{

q<-1

while
(q<=dim(t)[1])

{

png(paste(t[p,1],"-",t[q,1],'.',bin,'.heatmap.correlation.png',sep=''), width = 800, height = 800, units = "px", pointsize = 12)

heatmap.2 ((matcor[t[p,2]:t[p,3],t[q,2]:t[q,3]]), Colv=NA,Rowv=NA, dendrogram=c("none"),trace=c("none"), col=orgPal, ylab = t[p,1], xlab =t[q,1], labRow = NA, labCol = NA)

dev.off()

data_matrix_chr<-(data_matrix[t[p,2]:t[p,3],t[q,2]:t[q,3]])

matcor_chr<-matcor[t[p,2]:t[p,3],t[q,2]:t[q,3]]

rownames(matcor_chr) <- rownames(data_matrix_chr)

colnames(matcor_chr) <- colnames(data_matrix_chr)

write.table(matcor_chr,file = paste(t[p,1],"-",t[q,1],'.correl.txt',sep=''),sep = "\t")

q<-q+1

}

p<-p+1

}

setwd('../')







