
paste(OUTPUT,paste(method,a0,MAP.nontrivial[[j]],"RData",sep="."),sep=""))

load(paste(OUTPUT,paste(method,a0,MAP.nontrivial[[j]],"RData",sep="."),sep=""))

W <- res[[1]]
H <- res[[2]]
W1 <- W
H1 <- H
W.norm <- apply(W,2,function(x) x/sum(x))
for (i in 1:ncol(W)) {
  H1[i,] <- H1[i,]*colSums(W)[i]
  W1[,i] <- W1[,i]*rowSums(H)[i]
}
lambda <- res[[5]]
df <- get.df.solution(W1,H,lambda,tumor.type)
K <- length(table(df$class))

############# Signature plot
width <- 16
height <- ifelse(K==1,3,K*2)
pdf(file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"signature.pdf",sep="."),sep=""),width=width,height=height)
p <- plot.lego.observed.barplot(df,paste("Mutation Signatures in ",tumor.type,sep=""))
plot(p)
dev.off()

index.nonzero <- colSums(W) != 0
lambda <- unlist(res[[5]][length(res[[5]])])
lambda <- lambda/min(lambda)
lambda <- lambda[index.nonzero]
names(lambda) <- paste("W",seq(1:length(lambda)),sep="")
K <- sum(index.nonzero)
if (K != 1) {
  W0.tumor <- W[,index.nonzero]
  H0.tumor <- H[index.nonzero,]
  x <- W0.tumor
  y <- H0.tumor
  x.norm <- apply(x,2,function(x) x/sum(x))
  W1.tumor <- x.norm
  for (j in 1:K) y[j,] <- y[j,]*colSums(x)[j]
  H1.tumor <- y
  y.norm <- apply(y,2,function(x) x/sum(x))
  H2.tumor <- y.norm
} else {
  stop("No non-trivial solutations; All simulations converged to a trivial solution with one signature")
}

############# Reconstructing the activity of signatures
W.mid <- W1.tumor
H.mid <- H1.tumor
H.norm <- H2.tumor
if (length(grep("__",colnames(H.mid))) != 0) {
  hyper <- colnames(H.mid)[grep("__",colnames(H.mid))]
  H.hyper <- H.mid[,colnames(H.mid) %in% hyper]
  H.nonhyper <- H.mid[,!(colnames(H.mid) %in% hyper)]
  sample.hyper <- colnames(H.hyper)
  sample.hyper <- sapply(sample.hyper,function(x) strsplit(x,"__")[[1]][[1]])
  unique.hyper <- unique(sample.hyper)
  n.hyper <- length(unique.hyper)
  x.hyper <- array(0,dim=c(nrow(H.hyper),n.hyper))
  for (i in 1:n.hyper) {
    x.hyper[,i] <- rowSums(H.hyper[,sample.hyper %in% unique.hyper[i]])
  }
  colnames(x.hyper) <- unique.hyper
  rownames(x.hyper) <- rownames(H.mid)
  H.mid <- (cbind(H.nonhyper,x.hyper))
  H.norm <- apply(H.mid,2,function(x) x/sum(x))
}
W.norm <- apply(W.mid,2,function(x) x/sum(x))

##########################################################
############# W.norm = extracted signatures normalized to one
############# H.mid = activity of signatures across samples (expected mutations associated with signatures)
############# H.norm = normalized signature activity 
##########################################################
WH <- list(W.norm,H.mid,H.norm)
save(WH,file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"WH.RData",sep="."),sep=""))

############# Activity plot
p1 <- plot.activity.barplot(H.mid,H.norm,1.0,tumor.type)
pdf(file = paste(OUTPUT,paste(method,a0,"activity.barplot1",K,"pdf",sep="."),sep=""),width=15,height=12)
plot(p1)
dev.off()

main <- paste("Cosine similarity;",tumor.type,sep="")
pdf(file = paste(OUTPUT,paste(method,a0,"signature.comprison.sanger",K,"pdf",sep="."),sep=""),width=4,height=6)
s1 <- 1.5
s2 <- 2.0
par(mfrow=c(1,1))
par(mar=c(0,2,1,1))
x <- W.mid[1:96,]
colnames(x) <- paste("W",c(1:K),sep="")
plot.heatmap.2(plot.W.correlation(sanger,x),T,T,"ward.D",main)
dev.off()



lego96Percent<-apply(lego96,2,function(x){x/sum(x)})
lego96Count<-lego96
TCGABRCASignature<-W.mid
colnames(TCGABRCASignature)<-c("S1|Apobec","S2|Signature3","S3|Signature10","S4|DNA_mismatch_repair","S5|Age_related")
TCGABrca_Hmid<-H.mid
colnames(TCGABrca_Hmid)<-c("S1|Apobec","S2|Signature3","S3|Signature10","S4|DNA_mismatch_repair","S5|Age_related")
TCGABrca_Hnor<-H.norm
colnames(TCGABrca_Hnor)<-c("S1|Apobec","S2|Signature3","S3|Signature10","S4|DNA_mismatch_repair","S5|Age_related")
write.csv(TCGABrca_Hmid,"TCGA_EUROPE_BRCA_CGASignature_Hmid_Count.csv")

write.csv(TCGABrca_Hnor,"TCGA_EUROPE_BRCA_CGASignature_Hmid_Normalized.csv")



