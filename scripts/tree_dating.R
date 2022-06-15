#DATING TREES USING ASSOCIATED TIP DATES


#----FUNCTIONS----
n.days <- function(year, month) {
  # @param year: string (number) or integer for full year, e.g., '1997'
  # @param month: string (number) or integer for month
  
  start <- as.Date(paste(year, month, '01', sep='-'))
  if (month == '12') {
    # increment to start of next year
    end <- as.Date(paste(as.integer(year)+1, '01', '01', sep='-'))
  } else {
    end <- as.Date(paste(year, as.integer(month)+1, '01', sep='-'))
  }
  return(as.integer(difftime(end, start)))
}

mid.date <- function(lo, hi){
  mid <- as.integer(as.Date(lo,origin = "1970-01-01")) + as.integer(as.Date(hi,origin = "1970-01-01")-as.Date(lo,origin = "1970-01-01"))/2
  return (mid)
}

date.to.decimal <- function(dt) {
  return (as.double(dt) / 365.25 + 1970)
}

root.to.tip <- function(tre, vect, x){
  rtdtree <- tryCatch(
    {rtdtree <- rtt(tre, vect, ncpu=x)} ,
    error = function(c){
      message("Faulty rtt error, still running")
      return (NULL)
    }
  )
  return (rtdtree)
}

require(ape)

args <- commandArgs()

if (length(args) != 3)

dfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/Dates/edit",full.names=TRUE)
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/6Trees", full.names=TRUE)
afolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/5_1_final", full.names=TRUE)


branch.lengths <- data.frame()
genetic.dists <- data.frame()
big.df <- data.frame(stringsAsFactors = F)
setwd("~/PycharmProjects/hiv-evolution-master/8_Dated_Trees")
for (n in 1:length(tfolder)){
  tre <- read.tree(tfolder[n])
  dates <- read.csv(dfolder[n], header=FALSE,stringsAsFactors = F)
  names(dates) <- c('accno', 'date')
  
  filename <- strsplit(strsplit(tfolder[n], "/")[[1]][7], "_CR")[[1]][1]
  print(filename)
  
  # ROOTING THE TREES WITH RTT ----------------------------------
  tre <- multi2di(tre) 
  
  #rearranges the temp file to be in the same order as the tree
  index <- match(tre$tip.label, dates$accno)
  dates <- dates[index,]
  
  #loading a temporary data frame of mid dates
  
  temp <- NULL
  for (x in 1:nrow(dates)){
    dtx <- get.range(dates$date[x])
    temp <- rbind(temp, data.frame(accno=dates$accno[x],
                                   acc.date=dates$date[x], min.date=dtx[1], max.date=dtx[2],
                                   min.ord=date.to.decimal(dtx[1]), max.ord=date.to.decimal(dtx[2])))
   
  }
  sample.times <- date.to.decimal(mid.date(temp$min.date, temp$max.date))
  sample.times <- sample.times - 1970
  
  
  #ROOTING THE TREE
  rtdtree <- root.to.tip(tre, sample.times, 5)
  if (is.null(rtdtree)){
    print(filename)
    print("TREE PROBLEM")
  }
  
  #ROOTED (UNDATED) TREE ANALYSIS 
  # --------------------------------------------
  ntips <- Ntip(tre)
  
  # count the number of instances that first column (node) corresponds to a second column number which is <= n (meaning it is a tip)
  numtips <- tabulate(tre$edge[,1][which(tre$edge[,2] <= ntips)])
  
  #determines which nodes contain cherries (returns vector with their integer positions)
  is.cherry <- which(numtips == 2)
  
  # construct data frame where each row corresponds to a cherry
  m <- sapply(is.cherry, function(a) { #will input the numbers of nodes containing 2 tips?
    edge.idx <- tre$edge[,1]==a  # FINDS THE EDGES (row #) corresponding with the parent node ; ap: select rows in edge matrix with parent node index
    
    c(a,   # index of the parent node
      which(edge.idx),
      t(     # transpose so we can concatenate this vector to i
        tre$edge[edge.idx, 2]    # column of tip indices
      )
    )
  })
  df <- data.frame(node.index=m[1,], edge1=m[2,], edge2=m[3,], tip1=m[4,], tip2=m[5,])
  
  df$tip1.label <- tre$tip.label[df$tip1]
  df$tip2.label <- tre$tip.label[df$tip2]
  df$tip1.len <- tre$edge.length[df$edge1]
  df$tip2.len <- tre$edge.length[df$edge2]
  
  indels <- df[,c(6:9)]
  indels$total.length <- indels$tip1.len + indels$tip2.len
  
  #TERMINAL BRANCH LENGTHS PLOT 
  branch.lengths <- rbind(branch.lengths, data.frame(subtype=rep(filename,nrow(indels)), length=indels$total.length))
  
  #CHERRY GENETIC DISTANCES PLOT
  genetic.dists <- rbind(genetic.dists, data.frame(subtype=rep(filename,nrow(indels)), indels))
  
  #ROOT TO TIP DISTANCE PLOT
  sample.times <- sample.times + 1970
  lens <- node.depth.edgelength(rtdtree)[1:Ntip(rtdtree)]
  
  
  # STUFF USED FOR TEMPEST
  #rtdtree2 <- rtdtree
  #rtdtree2$tip.label <- paste(rtdtree2$tip.label, sample.times, sep="-")
  #setwd("~/PycharmProjects/hiv-evolution-master/7_Tempest")
  #write.tree(rtdtree2, file=paste0(filename,".tree"))
  
  big.df <- rbind(big.df,data.frame(subtype=rep(filename,Ntip(rtdtree)), dates=sample.times, lengths=lens))
  
  # DATING TREES USING LSD 
  # -----------------------------------
  write.tree(rtdtree, file="rtt2lsd.nwk")

  #The following is written in accordance with the LSD format

  write(nrow(temp), file='date_file.txt')

  for (i in 1:nrow(temp)) {
    if (temp$min.ord[i] == temp$max.ord[i]) {
       write(paste(temp$accno[i], temp$min.ord[i]), file='date_file.txt', append=TRUE)
     } else {
       write(paste(temp$accno[i], paste0("b(", temp$min.ord[i], ",", temp$max.ord[i], ")")),
             file='date_file.txt', append=TRUE)
     }
   }

  paths <- nodepath(rtdtree)
  distances <- c()
  for (p in 1:length(paths))  distances <- c(distances, sum(rtdtree$edge.length[paths[[p]]]))
  
  aLen <- length(read.FASTA(afolder[n])[1][[1]])
  
  #system(paste0('lsd -i rtt2lsd.nwk -d date_file.txt -o ', filename ,' -c -f 1000 -s ',aLen))
  
}


#genetic.dists <- cbind(genetic.dists, filtered.indels3[,7:21])
# BRANCH LENGTHS PLOT --------------------------------
#par(mar=c(6,7,2,2))
#branch.len2 <- split(branch.lengths$length, branch.lengths$subtype)
#boxplot(branch.len2, xlab="Subtype",cex.lab=1.3,las=1)
#par(las=3)
#mtext(side = 2, text = "Terminal Branch Lengths (Expected Substitutions)", line = 4, cex=1.3)

write.csv(genetic.dists, "~/vindels/Pipeline_2_Within/genetic-dists.csv")

#colnames(treedata) <- c("Clade","Estimate", "$x$-Intercept (Year)", "$R^2$")

#lxtable <- xtable(treedata)
#toLatex.xtable(lxtable)


#setwd("~/vindels/Indel_Analysis")
#write.csv(treedata,file="tree_data.csv")

# for (m in 1:7){
#   par(mar=c(4.5,6,1,1))
#   buffer <- 0.01*(range(new.df[m][[1]][,1])[2] - range(new.df[m][[1]][,1])[1])
#   newx <- min(new.df[m][[1]][,1])-buffer
#   plot(x=jitter(new.df[m][[1]][,1],amount=0.45), y=new.df[m][[1]][,2], col=alpha(cols[m],0.35), xlim=c(1980,2016),
#        ylim=c(0,0.30), cex=1,pch=20, 
#        xlab="Collection Date",ylab="Root-to-Tip Branch Length", cex.lab=1.1)
#   abline(lm(new.df[m][[1]][,2]~new.df[m][[1]][,1]),lwd=1.5)
#   text(newx+buffer*2,0.285,labels=names[m], cex=1.4)
#   par(xpd=NA)
#   text(newx-buffer*25,0.31,labels=paste0(letters[m],")"), cex=1.5)
#   par(xpd=F)
# }



