#install.packages(c("mixOmics","RVAideMemoire"))
library(mixOmics)
library(RVAideMemoire)

# Calculate whether measure is greater in control group after scaling
# Used to colour points in VIP plots
direction=function(dt,measure){
  dts = as.data.frame(scale(dt[,-1]))
  dts$Group = dt$Group
  res = median(dts[[measure]][dts$Group=="Control"],na.rm=TRUE) > median(dts[[measure]][dts$Group!="Control"],na.rm=TRUE)
  return(res)
}

or = rgb(247/255,158/255,84/255) # orange
bl = rgb(69/255,150/255,207/255) # blue

gr1 = rgb(116/255,161/255,77/255) # green
or1 = rgb(247/255,161/255,30/255) # orange

bl2 = rgb(3/255,40/255,210/255) # blue
or2 = rgb(68/255,186/255,82/255) # orange

# Read in data
dat = read.delim("Summary stats Control, patient and mouse 05-03-18.csv",sep=",",stringsAsFactors=FALSE)
rownames(dat) = dat$Subject
dat$Subject = NULL
rownames(dat)

# Transpose and ensure data are numeric
dat = data.frame(t(dat),stringsAsFactors=FALSE)
for(col in colnames(dat)[2:length(colnames(dat))]) dat[[col]]=as.numeric(as.character(dat[[col]]))

# Tidy measure names
colnames(dat) = gsub("MCI.Vol","MCIperVol",colnames(dat))
colnames(dat) = gsub("\\.","\n",colnames(dat))
colnames(dat) = gsub("\nmitos","",colnames(dat))
colnames(dat) = gsub("X\n","Percent",colnames(dat))
colnames(dat) = gsub("\nmito","",colnames(dat))
colnames(dat) = gsub("95\n\nCI","interval",colnames(dat))

# Pairwise correlation plot for examining relationship between measures
colours = c("black","red","blue")
names(colours) = c("Control", "Mito Disease", "Mouse")
pdat = dat[,-1]
pdf("pairwise.pdf",width=30,height=30,pointsize=17.5)
  pairs(pdat,col = colours[dat$Group],upper.panel = NULL, pch = 16, cex = 0.6)
dev.off()

# Consider dropping subset of measures
rejectroot = c("Mean","interval","Kurtosis","MCI","Vol","MCIperVol")
rejectmeasures = unlist(lapply(rejectroot,grep,colnames(dat)))

# Consider keeping subset of variables
keeproot = c("Nanotunnels\nper\n100","Volume\ndensity","Percent\nsimple","Percent\ncomplex","Percent\nsmall","Percent\nlarge","MCI\nMedian","Vol\nMedian")
keepmeasures = c(1,match(keeproot,colnames(dat)))

subsets = c("SUBSET1","SUBSET2","ALL")

# Consider two different comparisons: Control-Patient and Human(Control)-Mouse
keeps = list()
keeps[["Controls & Patients"]] = c("Control","Mito Disease")
keeps[["Controls & Mice"]] = c("Control","Mouse")

def.par = par(no.readonly = TRUE) # save default, for resetting...

for(subset in subsets){

  if(subset=="SUBSET2") datsub=dat[,keepmeasures]
  if(subset=="SUBSET1") datsub=dat[,-rejectmeasures]
  if(subset=="ALL") datsub=dat


  # Supervised learning with PLS-DA: find combinations of measures that best discriminate
  # between two groups.  PLS-DA also ranks original measures by their contribution to 
  # splitting categories (VIP score).  Number of components = 2, for ease of plotting.

  # Multi-page PDF report showing PLS-DA biplots and VIP scores
  pdf(paste("PLSDA",paste(subset,".pdf",sep=""),sep="_"),width=8.27*2.5,height=8.27,pointsize=18)
  layout(matrix(c(1,2,3),nrow=1,ncol=3,byrow=TRUE),widths=c(1,0.5,1))
  for (k in names(keeps)){
    keepcases = keeps[[k]]
    if ("Mouse" %in% keepcases){colcontrol=gr1; colother=or1}else{colcontrol=bl2; colother=or2}
    dt = datsub[datsub$Group%in%keepcases,]
    dt2 = dt[-6,]
    da = plsda(dt2[,-1], factor(dt2$Group), scale = TRUE, ncomp=2)
    op2=par(mar=c(4.5,4.5, 4, 2) + 0.1)
    plotIndiv(da,ellipse=TRUE,style="graphics",title=k,size.axis=1.5,size.xlabel=1.5,size.ylabel=1.5,col=c(colcontrol,colother),cex=1.5)
    par(op2)
    op3 = par(mar=c(4.5,7.0, 4, 2) + 0.1)
    tb = PLSDA.VIP(da)$tab
    if(subset=="ALL") {
       labs = gsub("\n"," ",rownames(tb))
       ulim=2
    }else{
       labs = gsub("\nper\n100","\nper 100",rownames(tb))
       ulim=1.5
    } 
    vcols = ifelse(sapply(rownames(tb),direction,dt = dt),colcontrol,colother)
    plot(tb$VIP,seq(length(tb$VIP),1),xlab="VIP",ylab="",axes=FALSE,type="n",xlim=c(0,ulim), cex.axis=1.5, cex.lab=1.5)
    abline(h = seq(length(tb$VIP),1),col="grey",lwd=2)
    abline(v=1,col="red",lty=2,lwd=2)
    points(tb$VIP,seq(length(tb$VIP),1),pch=15,cex=1.5,col=vcols[rownames(tb)])
    axis(1, cex.axis=1.5, cex.lab=1.5)
    axis(2, labels=labs,at = seq(length(tb$VIP),1),las=2, cex.axis=1.5, cex.lab=1.5)
    #legend("bottomright",legend=c("Higher in controls","Lower in controls"),col=c(colcontrol,colother),pch=15,bg="white")
    par(op3)
    plot(dt[[rownames(tb)[1]]],dt[[rownames(tb)[2]]],
      col=ifelse(dt$Group=="Control",colcontrol,colother),pch=16,cex=2,xlab=gsub("\n"," ",rownames(tb)[1]),ylab=gsub("\n"," ",rownames(tb)[2]),
      cex.axis=1.5, cex.lab=1.5)
    text(dt[[rownames(tb)[1]]],dt[[rownames(tb)[2]]],labels=rownames(dt),col=ifelse(dt$Group=="Control",colcontrol,colother),pos=4,cex=1.5)
    legend("topright",legend=keepcases,col=c(colcontrol,colother),pch=16,cex=1.5)
  }
  dev.off()

  # Unsupervised PCA also splits data nicely, however PCA results generally more difficult to 
  # interpret than PLS-DA for 2-way comparison.  Number of components = 2, for ease of plotting.

  # Multi-page PDF report showing PCA biplots and loading vectors
  pdf(paste("PCA",paste(subset,".pdf",sep=""),sep="_"),width=8.27*2,height=8.27)
  op = par(mfrow=c(1,2),mar=c(4,2, 4, 2) + 0.1)
  for (k in names(keeps)){
    keepcases = keeps[[k]]
    dt = datsub[datsub$Group%in%keepcases,]
    prc = pca(dt[,-1], scale = TRUE, center = TRUE, ncomp=2)
    biplot(prc,main=k,xlim=c(-0.6,0.6))
  }
  par(op)
  dev.off()
}