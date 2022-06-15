
# GENETIC DISTANCES PLOTS --------------

# qcquire the order of means 
new.df <- split(genetic.dists, genetic.dists$subtype)
means <- sapply(new.df, function (x) mean(x$GD) )
means <- means[order(means)]
ordered <- names(means)


gd.order <- genetic.dists
# this applies the properly ordered factor to the subtype column 
gd.order$subtype <- factor(gd.order$subtype, levels=ordered)
# this sorts the entire data frame by the properly ordered 'levels' in the subtype column
gd.order <- gd.order[order(gd.order$subtype),]


# A) Standard plot 
par(mar=c(5,5,1,1))
plot(gd.order$subtype, gd.order$GD, xlab="Group M Clade", ylab="Genetic Distance (Cherries)", cex.axis=1.2, cex.lab=1.5)

# perform a wilcoxon test on this 


# B) Only GDs under 0.05 
gd.05 <- genetic.dists[which(genetic.dists$GD <= 0.05),]
par(mar=c(5,5,1,1))
plot(gd.05$subtype, gd.05$GD, xlab="Group M Clade", ylab="Genetic Distance (Cherries)", cex.axis=1.2, cex.lab=1.5)


# C) Staggered density plot of the seven subtypes  
require(ggplot2)
require(ggridges)
gd.15 <- genetic.dists[which(genetic.dists$GD <= 0.15),]
ggplot(gd.15, aes(x=GD, y=subtype, group=subtype)) + 
  geom_density_ridges(colour="white", fill="blue3", scale=1, bandwidth=0.002) + 
  labs(x="Genetic Distance", y="Subtype") + 
  theme(axis.title.x=element_text(size=15,margin=margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=15,margin=margin(t = 15, r = 0, b = 0, l = 0)),
        #strip.text.x = element_text(size=16),
        axis.text=element_text(size=13))


big.df <- big.df[which(big.df$dates>1960),]
new.df <- split(big.df[,2:3],big.df[,1])

require(survey)
require(scales)
par(mar=c(5,5,4,3))

require(xtable)
require(RColorBrewer)
cols <- brewer.pal(7,"Dark2")
par(mfrow=c(4,2),las=1)
names <- c("AE", "AG","A1","B","C","D","F1")
letters <- c("a","b","c","d","e","f","g")
lims <- c()
treedata <- data.frame(stringsAsFactors = F)

for (m in 1:7){
  
  par(mar=c(4.5,6,1,1))
  xval <- new.df[m][[1]][,1]
  yval <- new.df[m][[1]][,2]
  
  #buffer <- 0.01*(range(xval)[2] - range(new.df[m][[1]][,1])[1])
  #newx <- min(new.df[m][[1]][,1])-buffer
  plot(x=jitter(xval,amount=0.45), y=yval, col=alpha(cols[m],0.35), xlim=c(1979,2016),
       ylim=c(0,0.30), cex=1,pch=20, 
       xlab="Collection Date",ylab="Root-to-Tip Branch Length", cex.lab=1.1)
  lreg <- lm(lengths~dates ,data=new.df[m][[1]])
  
  #compute the 95% conf ints for slope
  cislope <- c(strsplit(formatC(confint(lreg)[2],format="e",digits=2),"e")[[1]][1],
               strsplit(formatC(confint(lreg)[4],format="e",digits=2),"e")[[1]][1])
  rsqr <- signif(summary(lreg)$adj.r.squared,2)
  slope <- strsplit(formatC(lreg$coefficients[2][[1]],format="e",digits=2),"e")[[1]]
  
  
  #calculation of the xintercept
  xest <- svycontrast(lreg, quote(-`(Intercept)`/dates))
  xint <- coef(xest)[[1]]
  se <- SE(xest)[[1]]
  
  xlo <- format(round((xint-se),2),nsmall=2)
  xup <- format(round((xint+se),2),nsmall=2)
  treedata <- rbind(treedata, data.frame(Clade=names[m],
                                         Estimate= toString(sprintf("%s (%s,%s) %s",slope[1],cislope[1],cislope[2],slope[2])),
                                         xint=sprintf("%s (%s,%s)",format(round(xint,2),nsmall=2), xlo, xup),
                                         rsqr=rsqr,stringsAsFactors = F))
  
  
  
  abline(lreg,lwd=1.5)
  #print(summary(lreg))
  text(1979.7,0.285,labels=names[m], cex=1.6)
  eqn <- paste0("y = ", signif(lreg$coefficients[2][[1]],2), "x ", signif(lreg$coefficients[1][[1]],2))
  #text(1998,0.275,labels=eqn, cex=1.1)
  #rsqr <- paste0(txt, get(signif(summary(lreg)$adj.r.squared,2)))
  #text(2012,0.29,labels=expression(paste(R^2,"= ")), cex=1.1)
  #text(2015,0.286,labels=rsqr,cex=1.1)
  par(xpd=NA)
  text(1970,0.31,labels=paste0(letters[m],")"), cex=1.5)
  par(xpd=F)
  
  
}


#plot(x=jitter(dfAE$dates,amount=0.3), y=dfAE$lengths, col=alpha(cols[1],0.7), xlim=c(1982,2016),ylim=c(0.05,0.25), cex=0.8,pch=20,
#     xlab="Collection Date",ylab="Root-to-Tip Branch Length\n", main="Tree Branch Lengths vs. Collection Dates -- AE", cex.lab=1.4,cex.main=1.6)
#abline(lm(dfAE[,2]~dfAE[,1]),lwd=1.5)


#par(mar=c(6,6,4,4))
#plot(x=jitter(dfB$dates,amount=0.5), y=dfB$lengths, col=alpha(as.character(dfB$color),0.5), cex=0.8,pch=20,
#     xlab="Collection Date",ylab="Root-to-Tip Branch Length", main="Tree Branch Lengths vs. Collection Dates -- Subtype B",cex.lab=1.3,cex.main=1.6)
#abline(lm(dfB[,2]~dfB[,1]))

#par(mar=c(6,6,4,4))
#plot(x=jitter(dfC$dates,amount=0.5), y=dfC$lengths, col=alpha(cols[1],0.5),cex=1,pch=20, 
#     xlab="Collection Date",ylab="Root-to-Tip Branch Length", main="Tree Branch Lengths vs. Collection Dates -- Subtype C",cex.lab=1.3,cex.main=1.6)
#abline(lm(dfC[,2]~dfC[,1]))
