#Set file names/locations from command line
#Input file name required
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0){
    stop("The input file name must be supplied",call.=FALSE)
} else if (length(args)==1) {
    args[2] <- "paths"
}

#Function to plot CNVs
clustercnvsplot<-function(xdups,genelims,reglims=c(0,0),marfact=4,bin=1600){
    if (reglims[1]==0 & reglims[2]==0){  #If reglims not provided, set to the entire range of the file
        reglims<-range(as.numeric(colnames(xdups)))
    } else {
        temp<-xdups[,colnames(xdups) >reglims[1] & colnames(xdups)<reglims[2]]
        tempcnv<-temp[apply(temp,1,function(i) any(i %in% c(1,2,4,5))),,drop=FALSE]  #Only included individuals with copy number variation in the region specified
        xdups<-tempcnv
    }
    allcols<-c("#B8860B","#FFD12C","white","#88B3FF","#1185CC")
    presstates<-as.numeric(names(table(xdups)))
    cols<-allcols[min(presstates):max(presstates)] #assign colors for the states present in the file
    
    h <- hclust(dist(xdups))
    dd<-as.dendrogram(h,hang=-1)
    
    layout(matrix(c(1,2),nrow=1,byrow=TRUE), width = c(10,40))
    par(mar=c(0,1,0,7),cex=0.6)
    plot(dd,horiz=T,yaxt="n")
    
    par(mar=c(marfact,0,marfact,2))
    
    image(x = as.numeric(colnames(xdups))+as.numeric(bin)/2,
          y = c(1:nrow(xdups)),
          t(xdups[h$order,]),
          col = cols,
          axes=T,
          xlim=reglims,
          xlab="Position"
    )
    abline(v=genelims)
    legend("topleft",c("0","1","2","3","4"),fill=allcols,bty="n",title="Copy Number")
}

y2<-as.matrix(read.table(args[1],as.is=T,head=T,check.names=F))

#Plot CNVs
y2cnv<-y2[apply(y2,1,function(i) any(i %in% c(1,2,4,5))),]
pdf(sprintf("%s.plot.pdf",args[2]),height=10,width=8)
clustercnvsplot(y2cnv,bin=1600,marfact=3,genelims=c(144792018,144826716,144917256,144940496,145030455,145061904),reglims=c(144706830,145069066))
dev.off()
