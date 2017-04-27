library(msm)  #Library for truncated normal distribution

#Set file names/locations from command line
#Input file name required
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0){
    stop("The input file name must be supplied",call.=FALSE)
} else if (length(args)==1) {
    args[2] <- "paths"
}

#Function to calculate a bin-specific factor for how much above/below average coverage tends to be, using only individuals with copy number 2.
#Input to function is coverage (xcov), coverage means (mm), and current paths (yy), 
get_binfacts<-function(xcov,mm,yy){
    binfacts<-numeric()
    for (i in 1:ncol(xcov)){
        if (is.na(xcov[1,i])){
            binfacts[i]<-1
        } else {
            thisbin<-xcov[,i]/mm[[1]][,3]
            thisbin2<-thisbin[yy[,i]==3]
            binfacts[i]<-mean(thisbin2)
        }
    }
    return(binfacts)
}

#Function to update mu and sd using only bins with copy number 2
#Input to function is coverage (xcov), coverage means (mm), and current paths (yy), 
update_musd<-function(xcov,mm,yy){
    mus2<-mm[[1]]
    sds2<-mm[[2]]
    for (i in 1:nrow(xcov)){
        reg2<-which(yy[i,]==3 & !is.na(xcov[i,]))
        mus2[i,3]<-mean(xcov[i,reg2])
        sds2[i,3]<-sd(xcov[i,reg2])
    }
    mus2[,2]<-mus2[,3]*0.5
    mus2[,4]<-mus2[,3]*1.5
    mus2[,5]<-mus2[,3]*2
    sds2[,1]<-sds2[,3]*0.2
    sds2[,2]<-sds2[,3]*0.5
    sds2[,4]<-sds2[,3]*1.5
    sds2[,5]<-sds2[,3]*2
    return(list(mus2,sds2))
}

#Function to run viterbi
viterbi<-function(x,logtran,delta,xmsd,binfacts){
    ind<-nrow(x) #Number of individuals
    n<-ncol(x) #Length of data
    m <-nrow(logtran) #Number of states
    
    y<-list()
    
    emis<-list()
    logemis<-list()
    nu<-list()
    matrixnu<-list()
    
    xmus<-xmsd[[1]]
    xsds<-xmsd[[2]]
    
    uninform<-which(is.na(x[1,]))  #Indices of uninformative bins (denoted with NA; this assumes all individuals are uninformative for the same bins)
    
    for (i in 1:ind){
        #Set emission matrix for each individual, get emission probability of each state (cols) in each bin j (rows)
        emis[[i]]<-matrix(NA,nrow=n,ncol=m)
        for (j in setdiff(1:n,uninform)){
            emis[[i]][j,1]<-dtnorm(x[i,j],mean=xmus[i,1]*binfacts[j], sd=xsds[i,1], lower=0, upper=Inf)
            emis[[i]][j,2]<-dtnorm(x[i,j],mean=xmus[i,2]*binfacts[j],sd=xsds[i,2], lower=0, upper=Inf)
            emis[[i]][j,3]<-dtnorm(x[i,j],mean=xmus[i,3]*binfacts[j],sd=xsds[i,3], lower=0, upper=Inf)
            emis[[i]][j,4]<-dtnorm(x[i,j],mean=xmus[i,4]*binfacts[j],sd=xsds[i,4], lower=0, upper=Inf)
            emis[[i]][j,5]<-dtnorm(x[i,j],mean=xmus[i,5]*binfacts[j],sd=xsds[i,5], lower=0, upper=Inf)
        }
        emis[[i]][uninform,]<-c(1,1,1,1,1) #Set equal emission probabilities for all states in uninformative bins
        emis[[i]]<-emis[[i]]/rowSums(emis[[i]])
        logemis[[i]]<-log(emis[[i]])
        #Initialise
        nu[[i]]<-matrix(1,nrow=n,ncol=m)
        nu[[i]][1,]<-log(delta)+logemis[[i]][1,]
        
        #Fill traceback
        for (j in 2:n){
            matrixnu[[i]]<-matrix(nu[[i]][j-1,],nrow=m,ncol=m)
            nu[[i]][j,]<-apply(matrixnu[[i]]+logtran,2,max)+logemis[[i]][j,]
        }

        #Traceback
        y[[i]]<-rep(NA,n)
        y[[i]][n]<-which.max(nu[[i]][n,])
        for (k in seq(n-1,1,-1)){
            y[[i]][k]<-which.max(logtran[,y[[i]][k+1]]+nu[[i]][k,])
        }    
    }
    y2<-matrix(unlist(y), ncol = length(y[[1]]), byrow = TRUE)
    colnames(y2)<-colnames(x)
    rownames(y2)<-rownames(x)
    return(y2)
}

#Read binned coverage
allcov<-as.matrix(read.table(args[1],as.is=T,head=T,check.names=F))

#Set starting probabilities
delta<-c(0.01,0.01,1,0.01,0.01)
delta<-delta/sum(delta)

#Set transition matrix
tran<-matrix(c(0.9,0.00001,0.0001,0.00001,0.00001,0.00001,0.9,0.0001,0.00001,0.00001,0.001,0.001,1,0.001,0.001,0.00001,0.00001,0.0001,0.9,0.00001,0.00001,0.00001,0.0001,0.00001,0.9),nrow=5)
tran<-tran/rowSums(tran)
logtran<-log(tran)

#Set initial values for HMM
initmu2<-apply(allcov,1,function(z) mean(z,na.rm=T))    # Calculate mean for each individual
initsd2<-apply(allcov,1,function(z) sd(z,na.rm=T))      # Calculate sd for each individual

#Set states 0, 1, 3, 4 to have proportional means and sd; set state 0 sd to 0.2x
initmus<-initmu2 %*% matrix(c(0,0.5,1,1.5,2),1,5)
initsds<-initsd2 %*% matrix(c(0.2,0.5,1,1.5,2),1,5)
msd1<-list(initmus,initsds)

# Run viterbi iteratively, re-estimating bin factors and means and sd in between, until no further changes to paths
nochange<-FALSE
iteration<-0
y1<-matrix(nrow=nrow(allcov),ncol=ncol(allcov),data=3) #Initially set paths with no copy number variation
print("Running Viterbi algorithm")
while(!nochange){
    iteration<-iteration+1
    print(sprintf("Round %s",iteration))
    bb<-get_binfacts(allcov,msd1,y1) #Calculate bin factors
    y2<-viterbi(allcov,logtran,delta,msd1,binfacts=bb) #Run viterbi
    msd2<-update_musd(allcov,msd1,y2) # Calculate new means and sds.
    nochange<-all(y1 == y2) #Compare previous and new viterbi paths
    y1<-y2
    msd1<-msd2
}

write.table(y2,sprintf("%s.out.txt",args[2]),quote=F,sep="\t")

