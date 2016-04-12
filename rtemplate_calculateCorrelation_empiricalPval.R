library(ggplot2)
library(reshape2)
library(qvalue)

                                        # Script to compute correlation between tissues and TF enrichment across the 60 clusters
                                        # Script also calculates empirical p value after n rounds of random sampling and correlation calculation

args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)
d2 <- read.table(args[2], header=T)
n <- as.numeric(args[3])
## cluster motif overlap null_mean null_sd null_ks_d null_ks_pval enrich zscore pval qvalues significant

df <- d[,c("cluster","motif","enrich")]
d1 <- dcast(df, cluster~motif)

dnew <- merge(d1, d2)
dcorr <- as.data.frame(cor(dnew[c(2840:2870)], dnew[c(2:2839)]))
dcorr[is.na(dcorr)] <- 0
tissue <- rownames(dcorr)
dcorr <- data.frame(tissue, dcorr)

# to output table
## write.table(dcorr, file=args[3], quote=FALSE, sep='\t', row.names=FALSE)

ispositive <- function(x){
    x <- ifelse((x>0),1,0)
    return(x)
}

emp_pval_assign <- function(x){
    x <- 1
    return(x)
}

calc_emp_pval <- function(x){
    x <- x/(n+1)
    return(x)
}

d0 <- dnew[,2840:2870]
dscore <- data.frame(tissue, apply(dcorr[,-1], c(1,2), emp_pval_assign))

for (i in 1:n){
    dnew1 <- data.frame(dnew[,1:2839], d0[sample(nrow(d0)),])
    dcorr1 <- data.frame(tissue, as.data.frame(cor(dnew1[c(2840:2870)], dnew1[c(2:2839)])))
    dcorr1[is.na(dcorr1)] <- 0
    dscore[,-1] <- dscore[,-1] + apply((dcorr1[,-1] - dcorr[,-1]), c(1,2), ispositive)

    ## for (x in 1:a){
    ##     for (y in 2:b){
    ##         dscore[[x,y]] <- dscore[[x,y]] + isgreater(dcorr[[x,y]], dcorr$i[[x,y]])
    ##         ## print(x)
    ##         ## print(y)
    ##     }
    ## }
    ## ## print(dscore)
}
## dscore
## dscore[,-1] <- dscore[,-1]/n+1
dscore <- data.frame(tissue, apply(dscore[,-1], c(1,2), calc_emp_pval))
write.table(dscore, file=args[4], row.names=FALSE, quote=FALSE, sep='\t')

