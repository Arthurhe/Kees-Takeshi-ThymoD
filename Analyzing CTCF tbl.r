source("file:///F:/DATA/R/Kees/binding landscape/proB/CTCFanalysis_functionPack.R")
#############################
#ctcf tiling
#allctcfBED=read.table("file:///F:/DATA/R/Kees/Takeshi/loadingtable/all_peaks_annotation.bed",sep="\t")
#region_center_tiling(allctcfBED,"takeshi_all_ctcf_tiling10K.bed",bin_size = 10000)

#plot Bcl11b locus
target_region="Bcl11b"#IgH #mm9
chrom="chr12"
chromstart=103760000#112261789
chromend=110480000#116961527
##################################
CTCF_matrix_list_tot=c()
load(file="F:/DATA/R/Kees/Takeshi/CTCF_interaction_matrixs_KO_10K")#base0.05.Rdata
CTCF_matrix_list_tot[[1]]=CTCF_matrix_list
load(file="F:/DATA/R/Kees/Takeshi/CTCF_interaction_matrixs_WT_10K")
CTCF_matrix_list_tot[[2]]=CTCF_matrix_list
names(CTCF_matrix_list_tot)=c("KO","WT")
names(CTCF_matrix_list_tot$WT)=chr_list
names(CTCF_matrix_list_tot$KO)=chr_list
allctcfBED=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/all_CTCFpeak_directed.bed",sep="\t")
colnames(allctcfBED)=c("chr","start","end","id","strand","score")
rm(CTCF_matrix_list)

tagbedpe_WT=create_bedpe(allctcfBED[allctcfBED[,1]==chrom,],CTCF_matrix_list_tot$WT$chr12$PoisMatrix,threshold = 0.05)
tagbedpe_KO=create_bedpe(allctcfBED[allctcfBED[,1]==chrom,],CTCF_matrix_list_tot$KO$chr12$PoisMatrix,threshold = 0.05)

significant_interaction_WT=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/WT_CTCF_interaction.ctcf_idx_pair",sep="\t")
significant_interaction_KO=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/KO_CTCF_interaction.ctcf_idx_pair",sep="\t")
colnames(significant_interaction_WT)=c("idx1","idx2","significancy")
colnames(significant_interaction_KO)=c("idx1","idx2","significancy")
tagbedpe_WT=merge(tagbedpe_WT,significant_interaction_WT,by=c("idx1","idx2"))[,c(-1,-2)]
tagbedpe_KO=merge(tagbedpe_KO,significant_interaction_KO,by=c("idx1","idx2"))[,c(-1,-2)]

##############################
Exons_mm9=read.table("F:/DATA/R/Exons_mm9.sushibed",sep="\t",stringsAsFactors=F)
colnames(Exons_mm9)=c("chrom","start","stop","gene","score","strand","type")
Exons_mm9$strand[Exons_mm9$strand=="+"]=1
Exons_mm9$strand[Exons_mm9$strand=="-"]=-1
Exons_mm9$strand=as.numeric(Exons_mm9$strand)

#read PC1
WTpc1=read.table("file:///F:/DATA/R/Kees/Takeshi/wt.dn2.mm9.2r.40kR.80kSR.PC1.txt",sep="\t",skip=1)
WTpc1=WTpc1[,c(2,3,4,6)]
colnames(WTpc1)=c("chrom","start","stop","value")
KOpc1=read.table("file:///F:/DATA/R/Kees/Takeshi/84ko.dn2.mm9.3n.40kR.80kSR.PC1.txt",sep="\t",skip=1)
KOpc1=KOpc1[,c(2,3,4,6)]
colnames(KOpc1)=c("chrom","start","stop","value")
combindBedgraph=merge(WTpc1,KOpc1,by = c("chrom","start","stop"),all=T)
combindBedgraph=combindBedgraph[combindBedgraph[,1]==chrom&combindBedgraph[,2]>chromstart&combindBedgraph[,3]<chromend,]
pc1filp=(combindBedgraph[,4]>0)-(combindBedgraph[,5]>0)
pc1filp=c()
pc1filp[1]=combindBedgraph[88,2]
pc1filp[2]=combindBedgraph[138,3]

#precheck region
surrounding=100000
quest_region=Exons_mm9[Exons_mm9$gene==target_region,]
quest_region$start[1]=min(quest_region$start)-surrounding
quest_region$"stop"[1]=max(quest_region$"stop")+surrounding
quest_region=quest_region[1,]


#bedgraph reads
files1=list.files(path="F:/DATA/R/Kees/Takeshi/loadingtable/",full.names = T,pattern = ".*_Bcl11b.bedgraph")
names1=c("RNA1_KO","RNA2_KO","RNA3_KO","RNA1_WT","RNA2_WT","RNA3_WT")

combindBedgraph=NULL
for(i in 1:length(names1)){
  tmp=read.table(files1[i],sep="\t",stringsAsFactors=F)
  colnames(tmp)=c("chrom","start","stop","value")
  if(is.null(combindBedgraph)){combindBedgraph=tmp}
  else{combindBedgraph=merge(combindBedgraph,tmp,by = c("chrom","start","stop"),all=T)}
}
combindBedgraph[is.na(combindBedgraph)]=0
for(i in 1:length(names1)){
  tmp=combindBedgraph[,c(1,2,3,i+3)]
  colnames(tmp)[4]="value"
  assign(names1[i],tmp)
}

RNAreads_WT=cbind(RNA1_WT[,1:3],log2((RNA1_WT$value+RNA2_WT$value+RNA3_WT$value)/3+1))
RNAreads_KO=cbind(RNA1_KO[,1:3],log2((RNA1_KO$value+RNA2_KO$value+RNA3_KO$value)/3+1))
#methylation
methyl_WT=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/WT_Bcl11b_Methyllevel.bed",sep="\t")
methyl_WT[,4]=1-methyl_WT[,4]
methyl_KO=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/KO_Bcl11b_Methyllevel.bed",sep="\t")
methyl_KO[,4]=1-methyl_KO[,4]
#CTCF peaks
CTCFpeak_WT=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/WT_CTCF_160420_annotation.bed",sep="\t")
CTCFpeak_KO=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/lnc84KO_CTCF_160420_annotation.bed",sep="\t")
CTCFread_WT=read.table("file:///F:/DATA/R/Kees/Takeshi/loadingtable/WT_CTCF_Bcl11b_ucsc.bedgraph",sep="\t")
CTCFread_KO=read.table("file:///F:/DATA/R/Kees/Takeshi/loadingtable/lnc84KO_CTCF_Bcl11b_uscs.bedgraph",sep="\t")

threshold=0.01
#tagbedpeWT=create_bedpe(CTCFs_WT[CTCFs_WT[,1]==chrom,],CTCF_matrix_list_tot[[2]][[chrom]]$PoisMatrix,threshold = threshold)
#tagbedpeKO=create_bedpe(CTCFs_WT[CTCFs_WT[,1]==chrom,],CTCF_matrix_list_tot[[1]][[chrom]]$PoisMatrix,threshold = threshold)
tagbedpeWT=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/WTinteraction_Bcl11b.bedpe",sep="\t")
tagbedpeKO=read.table("F:/DATA/R/Kees/Takeshi/loadingtable/KOinteraction_Bcl11b.bedpe",sep="\t")
colnames(tagbedpeWT)=c("chrom1","start1","end1","chrom2","start2","end2","name","score") 
colnames(tagbedpeKO)=c("chrom1","start1","end1","chrom2","start2","end2","name","score") 

pdf(file="Takeshi_Bcl11b_locus_new.pdf",width=15,height=18)
yaxis="Normalized read depth"
CTCF_span=2000
#layout(matrix(c(1:3,4,6,1:3,5,7,8:10,11,13,8:10,12,14),5,4),heights = c(2.5,2.5,4,2.5,2.5))
#layout(matrix(c(rep(1,4),rep(12,4),rep(2,4),rep(13,4),rep(3,4),rep(14,4),4,4,5,5,15,15,16,16,6,6,7,7,17,17,18,18,8:11,19:22),6,8,byrow = T),heights = c(2.5,2.5,3.5,2.5,2.5,2.5))
layout(matrix(c(rep(1,6),rep(16,6),rep(2,6),rep(17,6),rep(3,6),rep(18,6),
                rep(4,3),rep(5,3),rep(19,3),rep(20,3),
                rep(6,3),rep(7,3),rep(21,3),rep(22,3),
                rep(8,3),rep(9,3),rep(23,3),rep(24,3),
                10:15,25:30),7,12,byrow = T),heights = c(1.5,2,3,2,2,2,2))
for(i in c(2,1)){
  currenttype=c("KO","WT")[i]
  #genebody
  pg = plotGenes(Exons_mm9[Exons_mm9$chrom==chrom & Exons_mm9$start>chromstart &Exons_mm9$stop< chromend,],chrom,chromstart,chromend ,
                 types = Exons_mm9[Exons_mm9$chrom==chrom & Exons_mm9$start>chromstart &Exons_mm9$stop< chromend,]$type, col="black", bentline = F,bheight=0.2,
                 labeltext=TRUE,maxrows=1,plotgenetype="box",fontsize = 0.6)
  labelgenome( chrom, chromstart,chromend,n=5,scale="Mb")
  #PC1
  plotBedgraph(get(paste(currenttype,"pc1",sep="")),chrom,chromstart,chromend,transparency=.50,color="gray0",range = c(-80,80))
  axis(side=2,las=2,tcl=.2)
  mtext("PC1 value",side=2,line=2.5,cex=0.6)
  labelgenome(chrom,chromstart,chromend,side=1,n=5,scale="Mb")
  title(main=paste(currenttype,"HiC PC1"))

  #interaction
  tagbedpe=get(paste("tagbedpe_",currenttype,sep=""))
  #colorlist=rep("gray20",nrow(tagbedpe))
  #colorlist[(tagbedpe[,"start1"]>pc1filp[1]&tagbedpe[,"start1"]<pc1filp[2])|(tagbedpe[,"start2"]>pc1filp[1]&tagbedpe[,"start2"]<pc1filp[2])]="firebrick4"
  colorlist=typetocol(tagbedpe$type)
  pbpe = plotBedpe(tagbedpe,chrom,chromstart,chromend,ylim=c(0,250),ylab="-log2 p-value",heights = -tagbedpe$significancy,plottype="loops",color=colorlist)
  axis(2)
  labelgenome(chrom,chromstart,chromend,side=1,n=5,scale="Mb")
  title(main=paste(currenttype,"significant Interaction"))
  rm(tagbedpe)
  legend("topleft",inset=0,legend=c("convergent","concurrent","divergent","unknown"),bty = "n",horiz =T,
         fill=c("firebrick4","dodgerblue4","gray0","gray"),border=c("firebrick4","dodgerblue4","gray0","gray"),text.font=2,cex=1.0)
  #CTCF
  plotRangesOver(get(paste("CTCFpeak_",currenttype,sep="")),chrom,chromstart,chromend,col="black",height=2,ybottom=0.1)
  #H33 zoom in #tagbedpe[tagbedpe$start1>85300000 &tagbedpe$end1<85500000,]
  zoomregion1 = c(108306419,108402000)
  zoomregion2 = c(109130000,109260000)
  zoomsregion(zoomregion1,extend=c(0.01,0.11),wideextend=0.11,offsets=c(0,0.560))
  zoomsregion(zoomregion2,extend=c(0.01,0.11),wideextend=0.11,offsets=c(0.560,0))
  ##the above part need to be tailed for each gene
  #CTCF
  plotBedgraph(get(paste("CTCFread_",currenttype,sep="")),chrom,chromstart=zoomregion1[1],transparency=.50,chromend=zoomregion1[2],color="firebrick2",range = c(-0.6,8))
  plotArrowheadOver(CTCFpeak_WT,chrom,chromstart,chromend,col="black",ybottom=-0.4,arr.length = 0.1,angle=90)
  labelgenome(chrom,chromstart=zoomregion1[1],chromend=zoomregion1[2],n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
  zoombox(passthrough=T)
  mtext(yaxis,side=2,line=2.5,cex=0.6)
  title(main=paste(currenttype,"CTCF reads in\nBcl11b enhancer region"))
  axis(side=2,las=2,tcl=.2)
  plotBedgraph(get(paste("CTCFread_",currenttype,sep="")),chrom,chromstart=zoomregion2[1],transparency=.50,chromend=zoomregion2[2],color="firebrick2",range = c(-0.6,8))
  plotArrowheadOver(CTCFpeak_WT,chrom,chromstart,chromend,col="black",ybottom=-0.4,arr.length = 0.1,angle=90)
  labelgenome(chrom,chromstart=zoomregion2[1],chromend=zoomregion2[2],n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
  zoombox(passthrough=T)
  mtext(yaxis,side=2,line=2.5,cex=0.6)
  title(main=paste(currenttype,"CTCF reads in\nBcl11b gene region"))
  axis(side=2,las=2,tcl=.2)
  #RNA
  plotBedgraph(get(paste("RNAreads_",currenttype,sep="")),chrom,chromstart=zoomregion1[1],transparency=.50,chromend=zoomregion1[2],color="gold2",range = c(-0.6,8))
  plotRangesOver(get(paste("CTCFpeak_",currenttype,sep="")),chrom,chromstart,chromend,col="black",height=1.6,ybottom=0)
  plotArrowheadOver(CTCFpeak_WT,chrom,chromstart,chromend,col="black",ybottom=-0.4,arr.length = 0.1,angle=90)
  labelgenome(chrom,chromstart=zoomregion1[1],chromend=zoomregion1[2],n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
  zoombox(passthrough=T)
  mtext(yaxis,side=2,line=2.5,cex=0.6)
  title(main=paste(currenttype,"RNA reads in\nBcl11b enhancer region"))
  axis(side=2,las=2,tcl=.2)
  plotBedgraph(get(paste("RNAreads_",currenttype,sep="")),chrom,chromstart=zoomregion2[1],transparency=.50,chromend=zoomregion2[2],color="gold2",range = c(-0.6,8))
  plotRangesOver(get(paste("CTCFpeak_",currenttype,sep="")),chrom,chromstart,chromend,col="black",height=1.6,ybottom=0)
  plotArrowheadOver(CTCFpeak_WT,chrom,chromstart,chromend,col="black",ybottom=-0.4,arr.length = 0.1,angle=90)
  labelgenome(chrom,chromstart=zoomregion2[1],chromend=zoomregion2[2],n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
  zoombox(passthrough=T)
  mtext(yaxis,side=2,line=2.5,cex=0.6)
  title(main=paste(currenttype,"RNA reads in\nBcl11b gene region"))
  axis(side=2,las=2,tcl=.2)
  #Methyl
  plotBedgraph(get(paste("methyl_",currenttype,sep="")),chrom,chromstart=zoomregion1[1],transparency=.50,chromend=zoomregion1[2],color="dodgerblue2",range = c(-0.075,1))
  plotRangesOver(get(paste("CTCFpeak_",currenttype,sep="")),chrom,chromstart,chromend,col="black",height=0.2,ybottom=0)
  plotArrowheadOver(CTCFpeak_WT,chrom,chromstart,chromend,col="black",ybottom=-0.05,arr.length = 0.1,angle=90)
  zoombox(passthrough=F)
  mtext("Hypomethylation level",side=2,line=2.5,cex=0.6)
  title(main=paste(currenttype,"Hypomethylation levelin \nBcl11b enhancer region"))
  axis(side=2,las=2,tcl=.2)
  labelgenome(chrom,chromstart=zoomregion2[1],chromend=zoomregion2[2],n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
  #setup zoom
  zoomregion3 = c(108324501-CTCF_span,108324773+CTCF_span)
  zoomregion4 = c(108356586-CTCF_span,108356858+CTCF_span)
  zoomregion5 = c(108387164-CTCF_span,108387436+CTCF_span)
  zoomsregion(zoomregion3,extend=c(0.05,0.3),wideextend=0.2,offsets=c(0,0.853))
  zoomsregion(zoomregion4,extend=c(0.05,0.3),wideextend=0.2,offsets=c(0.426,0.426))
  zoomsregion(zoomregion5,extend=c(0.05,0.3),wideextend=0.2,offsets=c(0.853,0))
  
  plotBedgraph(get(paste("methyl_",currenttype,sep="")),chrom,chromstart=zoomregion2[1],transparency=.50,chromend=zoomregion2[2],color="dodgerblue2",range = c(-0.075,1))
  plotRangesOver(get(paste("CTCFpeak_",currenttype,sep="")),chrom,chromstart,chromend,col="black",height=0.2,ybottom=0)
  plotArrowheadOver(CTCFpeak_WT,chrom,chromstart,chromend,col="black",ybottom=-0.05,arr.length = 0.1,angle=90)
  zoombox(passthrough=F)
  mtext("Hypomethylation level",side=2,line=2.5,cex=0.6)
  title(main=paste(currenttype,"Hypomethylation levelin \nBcl11b enhancer region"))
  axis(side=2,las=2,tcl=.2)
  labelgenome(chrom,chromstart=zoomregion2[1],chromend=zoomregion2[2],n=4,scale="Kb",edgeblankfraction=0.2,cex.axis=.75)
  #setup zoom
  zoomregion6 = c(109215286-CTCF_span,109215558+CTCF_span)
  zoomregion7 = c(109242248-CTCF_span,109242520+CTCF_span)
  zoomregion8 = c(109244438-CTCF_span,109244710+CTCF_span)
  zoomsregion(zoomregion6,extend=c(0.05,0.3),wideextend=0.2,offsets=c(0,0.853))
  zoomsregion(zoomregion7,extend=c(0.05,0.3),wideextend=0.2,offsets=c(0.426,0.426))
  zoomsregion(zoomregion8,extend=c(0.05,0.3),wideextend=0.2,offsets=c(0.853,0))
  
  #H33 zoom in #tagbedpe[tagbedpe$start1>85300000 &tagbedpe$end1<85500000,]
  for(i in 3:8){
    tagzoom=get(paste("zoomregion",i,sep=""))
    plotBedgraph(get(paste("methyl_",currenttype,sep="")),chrom,chromstart=tagzoom[1],transparency=.50,chromend=tagzoom[2],color="dodgerblue2",range = c(-0.075,1))
    plotRangesOver(get(paste("CTCFpeak_",currenttype,sep="")),chrom,chromstart,chromend,col="black",height=0.2,ybottom=0)
    plotArrowheadOver(CTCFpeak_WT,chrom,chromstart=tagzoom[1],chromend=tagzoom[2],col="black",ybottom=-0.05,arr.length = 0.1,angle=90)
    zoombox(passthrough=F)
    mtext("Hypomethylation level",side=2,line=2.5,cex=0.6);axis(side=2,las=2,tcl=.2)
    #labelgenome(NULL,chromstart=tagzoom[1],chromend=tagzoom[2],n=1,scale="Kb",edgeblankfraction=0.2,cex.axis=.75,chromcex=0.75, scalecex=0.75 )
    #abline(h=0.5,lty=2)
  }
}
dev.off()



