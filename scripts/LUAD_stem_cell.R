dir.create('scripts')
dir.create('results')
dir.create('fig')
rm(list = ls())
options(stringsAsFactors = F)

plot_cor_point=function(x,y,method='Pearson',top_col='#D55E00',right_col='#009E73'
                        ,ylab='y expression',xlab='x expression',title=NULL
                        ,marginal.type=c("histogram", "boxplot", "density", "violin", "densigram")[1]){
  library(ggplot2)
  library(ggpubr)
  library(ggExtra)
  corT=cor.test(x,y,method="pearson")
  cor=corT$estimate
  pValue=corT$p.value
  df=data.frame(x,y)
  p=ggplot(df, aes(x, y)) + 
    xlab(xlab)+ylab(ylab)+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'pearson', aes(x =x, y =y))+labs(title = title)
  p=ggMarginal(p, type = marginal.type, xparams = list(fill = top_col),yparams = list(fill = right_col))
  return(p)
}

wb_boxplot <- function(dat = data,
                       groups = groups,
                       xlab = '',
                       ylab = '',
                       xangle = 90,
                       title = 'Groups',
                       col = mycolor) {
  tmp.dat <- data.frame()
  for (ge in rownames(dat)) {
    print(ge)
    tmp <- data.frame(Samples = colnames(dat),
                      Genes = ge,
                      Values = as.numeric(dat[ge, ]),
                      Groups = groups)
    tmp.dat <- rbind(tmp.dat, tmp)
  }
  
  
  print(head(tmp.dat))
  library(ggpubr)
  if (length(unique(groups)) > 2) {
    tmp_plot <- ggplot(tmp.dat, 
                       aes(x=Genes, y=Values, fill=Groups)) +
      geom_boxplot(notch = F) +  
      stat_compare_means(method = "anova", label = "p.signif") +
      scale_fill_manual(values = col) + theme_test() +
      theme(axis.text.x = element_text(angle=xangle, 
                                       hjust = 0.95,
                                       vjust = 0.95)) +
      xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(title = title))
    return(tmp_plot)
  } else {
    tmp_plot <- ggplot(tmp.dat, 
                       aes(x=Genes, y=Values, fill=Groups)) +
      geom_boxplot(notch = F) +  
      stat_compare_means(method = "t.test", label = "p.signif") +
      scale_fill_manual(values = col) + theme_test() +
      theme(axis.text.x = element_text(angle=xangle, 
                                       hjust = 0.95,
                                       vjust = 0.95)) +
      xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(title = title))
    return(tmp_plot)
  }
}
ggplotKMCox=function(dat,title='Groups',labs=NULL,add_text=NULL, col = mycolor){
  library(ggplot2)
  colnames(dat)=c('time','status','groups')
  #sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #print((sdf))
  #summary(sdf)
  #p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  surv=survminer::ggsurvplot(sf, data = dat, palette = col, #jco palette 
                             pval = TRUE, surv.median.line='hv',
                             linetype = "strata"
                             ,conf.int = T
                             # ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             legend.title = title
                             ,legend.labs = labs
  )
  p1=surv$plot+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                ,axis.text.x=element_blank()
                                ,axis.title.x=element_blank()
                                ,plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches")
                                #,axis.title.y=element_blank()
                                ,legend.position=c(1,1), legend.justification=c(1,1)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  #p1=p1+text()
  #tms=data.frame(Group=tms.gp,value=tms.tps,Attribute=rep(data_m[1,1],length(tms.gp))
  #               ,ymax=rep(max(ylim),length(tms.gp)))
  #p4=p4+geom_text(data=tms,aes(x=Group, y=ymax, label=value),color="black")
  if(!is.null(add_text)){
    text.tb=surv$data.survplot[1,]
    text.tb[1,1]=0
    text.tb[1,5]=0
    text.tb$Text=add_text
    p1=p1+geom_text(data=text.tb,aes(x=time, y=surv, label=Text),color="black",hjust =0)
  }
  
  p2=surv$table+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                 #,axis.text.x=element_blank()
                                 #,axis.title.x=element_blank()
                                 #,axis.title.y=element_blank()
                                 ,plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches")
                                 ,plot.title=element_blank()
                                 ,legend.position=c(1,1), legend.justification=c(1,1)
                                 #,legend.background = element_rect(fill = NA, colour = NA)
                                 ,legend.title = element_text(family="Times",face="plain")
                                 ,legend.text = element_text(family="Times",face="plain"))
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  return(g2)
}
mg_volcano_wb <- function(logfc,
                          pvalue,
                          symbol=NULL,
                          cutFC=1,
                          cutPvalue=0.05
                          ,showText=NULL
                          ,colors=c(mycolor[2],'black',mycolor[1])
                          ,xlim=NULL,ylim=NULL
                          ,legend.pos='right'
                          ,ylab='-log10(FDR)',
                          leg='State',
                          xlab='log2(FoldChange)'){
  library(ggplot2)
  pos=c(0,0)
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else{
    pos='right'
  }
  cange=rep('None',length(logfc))
  cange[which(logfc>cutFC&pvalue<cutPvalue)]='Up'
  cange[which(logfc< -cutFC&pvalue<cutPvalue)]='Down'
  if(is.null(symbol)){
    symbol=rep('',length(logfc))
    showText=NULL
  }
  vo.input=data.frame(logFC=logfc,FDR=pvalue,change=cange,SYMBOL=symbol)
  #print(head(vo.input))
  p1 <- ggplot(data = vo.input, 
               aes(x = logFC, 
                   y = -log10(FDR)))
  p1=p1+geom_point(alpha=0.85, size=2, aes(color=change))  
  p1=p1+scale_color_manual(values=colors,limits = c("Down",'None', "Up"),name=leg) 
  p1=p1+geom_vline(xintercept=c(-cutFC,cutFC),lty=6,col="black",lwd=0.8)  
  p1=p1+geom_hline(yintercept = -log10(cutPvalue),lty=6,col="black",lwd=0.8)  
  p1=p1+ylab(ylab)+xlab(xlab)
  p1=p1+theme_bw()
  p1=p1+theme(
    axis.text.y=element_text(family="Times",face="plain"), #
    axis.title.y=element_text(family="Times",face="plain"), #
    legend.text=element_text(face="plain", family="Times", colour="black"  #
    ),
    legend.title=element_text(face="plain", family="Times", colour="black" #
    ),
    legend.justification=pos, legend.position=pos
    ,legend.background = element_rect(fill = NA, colour = NA)
  )
  if(is.null(showText)|is.null(symbol)){
    showText=c()
  }
  
  if(length(showText)>0){
    for_label <-vo.input[match(intersect(showText,vo.input$SYMBOL),vo.input$SYMBOL),]
    p1=p1+geom_point(size = 3, shape = 1, data = for_label)+ggrepel::geom_label_repel(
      aes(label = SYMBOL),
      data = for_label,
      color="black"
    )
  }
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(!is.null(xlim)){
    p1=p1+xlim(xlim)
  }
  
  return(p1)
}
mg_violin <- function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',legend.pos=NULL,melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,mycolor = mycolor){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=1)+scale_color_manual(values = mycolor)
  if(ct<=4){
    p1=p1+ggsci::scale_fill_d3()
  }else if(ct<=10){
    p1=p1+ggsci::scale_fill_jco(name=leg.title)
  }else if(ct<=20){
    p1=p1+ggsci::scale_fill_jco(palette = "category20",name=leg.title)
  }else if(ct<=30){
    cbPalette=c(ggsci::pal_jco("nrc", alpha = 0.6)(10),ggsci::pal_jco("category20", alpha = 0.6)(20))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }else if(ct<=38){
    cbPalette=c(ggsci::pal_d3()(10)
                ,ggsci::pal_jco("nrc", alpha = 0.6)(10)
                ,ggsci::pal_jco("category20", alpha = 0.6)(20)
                ,ggsci::pal_jco("default", alpha = 0.6)(8))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }
  
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif", step_increase = 0.0)
  }
  return(p1)
}
ggplotTimeROC=function(time,status,score,mks=c(1,3,5), col = mycolor){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    #roc.tm=ROC.DSST$times[which(ROC.DSST$times>0)]
    
    #p.dat=rbind()
    #for(i in which(ROC.DSST$times>0)){
    #los=lowess(ROC.DSST$FP[,i], y=ROC.DSST$TP[,i], f = 1/3, iter = 100)
    #los$x=c(0,los$x,1)
    #los$y=c(0,los$y,1)
    # p.dat=rbind(p.dat,data.frame(los$x, y=los$y,rep(lbs[i],length(los$y)),stringsAsFactors = F))
    #}
    
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+
      theme_bw()+
      xlab('False positive fraction')+
      ylab('True positive fraction') +
      scale_color_manual(values = col)
    #p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction') 
    
    p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_text(family="Times",face="plain")
                ,axis.title.x=element_text(family="Times",face="plain"),axis.title.y=element_text(family="Times",face="plain")
                ,plot.title=element_blank()
                ,plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches")
                ,legend.position=c(1,0)
                ,legend.justification=c(1,0)
                ,legend.background = element_rect(fill = NA, colour = NA)
                ,legend.title = element_text(family="Times",face="plain")
                ,legend.text = element_text(family="Times",face="plain"))
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}
mg_PlotMutiBoxplot=function(data,group,group_cols='jco'
                            ,test_method=c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[1]
                            ,order=NULL,size=1,fill=F,outlier.shape=NA,yscale=c('none','log2','log10')[1]
                            ,xangle=45,ylab='Value',xlab='',box_span=0.7
                            ,orientation = c("vertical", "horizontal", "reverse")[1]
                            ,legend.pos=NULL,melt=F,ylim=NULL,binwidth=0.05
                            ,add=c("none", "dotplot", "jitter", "boxplot", "point", "mean"
                                   , "mean_se", "mean_sd", "mean_ci", "mean_range", "median"
                                   , "median_iqr", "median_mad", "median_range")[3]){
  paired=FALSE
  if(test_method=='paired_t.test'|test_method=='paired_wilcox.test'){
    test_method=gsub('paired_','',test_method)
    paired=TRUE
  }
  print(class(data))
  if(add=='jitter'){
    fill=F
  }
  
  library(ggplot2)
  if(!melt){
    #print(class(data))
    if(class(data)=='numeric'|class(data)=='integer'){
      data=as.numeric(data)
      vd1.sbs=data.frame(group,rep('Tag',length(group)),data)
      #print(vd1.sbs)
    }else{
      data=as.data.frame(data)
      data$ID=group
      vd1.sbs <- reshape2::melt(data, id.vars=c("ID"))
    }
    colnames(vd1.sbs)=c('category','type','Score')
    Data=vd1.sbs
  }else{
    vd1.sbs=data
    colnames(vd1.sbs)=c('category','type','Score')
    Data=vd1.sbs
  }
  
  #vd1.sbs[,2]=paste0('C',as.numeric(as.character(vd1.sbs[,2])))
  if(is.null(order)){
    order=unique(vd1.sbs[,2])
  }
  
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos=c(0,0)
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='top'){
    pos='top'
  }else if(legend.pos=='buttom'){
    pos='buttom'
  }else{
    pos='right'
  }
  print(pos)
  if(fill){
    p <- ggpubr::ggboxplot(vd1.sbs, x="type", y="Score", fill = "category", yscale = yscale
                           ,palette = group_cols,width = box_span,size = size,order = order,outlier.shape=outlier.shape
                           ,orientation=orientation,add=add,add.params = list(binwidth=binwidth)
                           ,short.panel.labs = T)#
  }else{
    p <- ggpubr::ggboxplot(vd1.sbs, x="type", y="Score", color = "category", yscale = yscale
                           ,palette = group_cols,width = box_span,size = size,order = order,outlier.shape=outlier.shape
                           ,orientation=orientation,add=add,add.params = list(binwidth=binwidth)
                           ,short.panel.labs = T)#
  }
  
  p=p+ggpubr::stat_compare_means(aes(group=category), label = "p.signif", method = test_method,paired=paired
                                 #,label.y = max(vd1.sbs[,1])
  )
  #p=p+ylab(ylab)+xlab(xlab)
  #p=p+theme(axis.text.x = element_text(angle = xangle, hjust = 1))
  p=p+theme_bw()+theme(axis.text.x=tx, #
                       axis.text.y=element_text(family="Times",face="plain"), #
                       axis.title.y=element_text(family="Times",face="plain"), #
                       #panel.border = element_blank(),axis.line = element_line(colour = "black"), #
                       legend.text=element_text(face="plain", family="Times", colour="black"  #
                       ),
                       legend.title=element_text(face="plain", family="Times", colour="black" #
                       ),
                       legend.justification=pos, legend.position=pos
                       ,legend.background = element_rect(fill = NA, colour = NA)
                       #,panel.grid.major = element_blank(),   #
                       #panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab) #
  if(!is.null(ylim)){
    p=p+ylim(ylim)
  }
  return(p)
}

plotMutiBar=function(dat,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) + theme_bw()
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_jco()+theme(legend.position = "bottom")
  pg=pg+ggsci::scale_fill_jco()+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") + theme_bw()+xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(anova p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
  if(is.null(lambda)){
    lmda=cv_fit$lambda.min
  }else{
    lmda=lambda
  }
  fit.coef=fit$beta[(apply(fit$beta,1,function(x){
    return(sum(x!=0))
  })>0),]
  
  fit.coef=as.matrix(fit.coef)
  colnames(fit.coef)=fit$lambda
  #fit$lambda==cv_fit$lambda
  library(ggplot2)
  dat=data.table::melt(t(as.matrix(fit.coef)))
  dat_z=dat[which(dat$value==0),]
  dat=dat[which(dat$value!=0),]
  dat.sv=rbind()
  for (u in unique(dat_z[,2])) {
    t.z=dat_z[which(dat_z[,2]==u),1]
    t.zx=max(t.z)
    dat.sv=rbind(dat.sv,c(t.zx,u,0))
    t.zn=min(t.z)
    if(t.zx!=t.zn){
      dat.sv=rbind(dat.sv,c(t.zn,u,0))
    }
  }
  colnames(dat.sv)=colnames(dat_z)
  #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
  dat=crbind2DataFrame(rbind(dat,dat.sv))
  mn=min(-log(dat$Var1))
  mx=max(-log(dat$Var1))
  if(show_text){
    mx=(mx-mn)*0.1+mx
  }
  p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
  p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
  if(show_text){
    fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
    for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
    p=p+ggrepel::geom_label_repel(
      aes(label = Var2,color=Var2),
      data = for_label,hjust = 0
    )
  }
  p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
  p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
  tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                 ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
  p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
    geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
    geom_point(aes(colour=col))
  p1=p1+theme_bw()+theme(legend.position = "none")
  gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                        #,align = "hv"
                        ,labels = figLabels)
  return(gal)
}

cor_point <- function(x,
                      y,
                      method='Pearson',
                      top_col='red',
                      right_col='blue',
                      ylab='y expression',
                      xlab='x expression'){
  #x=rnorm(200)
  #y=rnorm(200)
  library(ggplot2)
  data=data.frame(x,y)  
  colnames(data)=c('wt','mpg')
  til=''
  if(method=='Pearson'){
    cr=cor.test(x,y)
    p=cr$p.value
    r=cr$estimate
    til=paste0('Pearson\'s correlation\nR=',round(r,3),'\nP=',signif(p,3))
  }else{
    cr=cor.test(x,y,method = "spearman")
    p=cr$p.value
    r=cr$estimate
    til=paste0('spearman correlation\nR=',round(r,3),'\nP=',signif(p,3))
  }
  
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(                              
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
      ,plot.margin=unit(c(0.1, 0.1, 0, 0), "inches")
    )
  empty=empty+geom_text(aes(x=1, y=1, label=til),color="black")
  
  plot_top <- ggplot(data, aes(wt, fill=top_col)) + 
    geom_density(alpha=.5,fill=top_col) +theme_bw()+ 
    theme(legend.position = "none",                           
          #plot.background = element_blank(), 
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          #panel.background = element_blank(),
          axis.title.x = element_blank(),
          #axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks.x = element_blank()
          ,plot.margin=unit(c(0.1, 0, 0, 0.1), "inches")
    )
  
  plot_right <- ggplot(data, aes(mpg, fill=right_col)) + 
    geom_density(alpha=.5,fill=right_col) +coord_flip()+theme_bw()+ 
    theme(legend.position = "none",                           
          #plot.background = element_blank(), 
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          #panel.background = element_blank(),
          #axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
          ,plot.margin=unit(c(0.01, 0.1, 0.1, 0), "inches")
    )
  #scale_fill_manual(values = c("orange", "purple")) + 
  
  p1=ggplot(data=data, aes(x=wt, y=mpg))+geom_point()+stat_smooth(method="lm")
  p1=p1+theme_bw()
  p1=p1+theme(axis.text.x=element_text(family="Times",face="plain"), #
              axis.text.y=element_text(family="Times",face="plain"), #
              axis.title.y=element_text(family="Times",face="plain"), #
              #panel.border = element_blank(),
              axis.line = element_line(colour = "black"), #
              legend.text=element_text(face="plain", family="Times", colour="black"),  #
              legend.title=element_text(face="plain", family="Times", colour="black"), #
              plot.margin=unit(c(0.01, 0.01, 0.1, 0.1), "inches")
              #,panel.grid.major = element_blank(),   #
              #panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab)
  
  pg1=ggpubr::ggarrange(plot_top,p1, ncol = 1, nrow = 2,heights = c(0.3,1),align = "v")
  pg2=ggpubr::ggarrange(empty,plot_right, ncol = 1, nrow = 2,heights = c(0.3,1),align = "v")
  
  pgal=ggpubr::ggarrange(pg1,pg2, ncol = 2, nrow = 1,widths = c(1,0.3),align = "h")
  return(pgal)
}

plot_GSEA_By_nodes_wb <- function(parseGSEAResult,
                                  indexs=c(1,2,3),
                                  color = mycolor,
                                  TermNames=NULL,
                                  left=NULL,
                                  right=NULL){
  #parseGSEAResult=gsea.result.kegg.result
  library(ggplot2)
  if(is.null(parseGSEAResult$TEMPLATE)){
    if(is.null(left)){
      left='RankTop'
    }
    if(is.null(right)){
      right='RankBottom'
    }
  }
  if(is.null(left)){
    left=parseGSEAResult$TEMPLATE[1]
  }
  if(is.null(right)){
    right=parseGSEAResult$TEMPLATE[2]
  }
  if(!is.null(TermNames)){
    inds=match(TermNames,parseGSEAResult$EnrichTable[,1])
    inds=inds[!is.na(inds)]
    if(length(inds)==0){
      print(paste0(TermNames,' Not Found!'))
      return(NA)
    }
  }else{
    inds=indexs
    if(max(inds)>nrow(parseGSEAResult$EnrichTable)){
      print(paste0(inds,' out range!'))
      return(NA)
    }
  }
  #parseGSEAResult=GSE17705_GSEA
  g.rnk=parseGSEAResult$Rank
  
  all.info=rbind()
  all.dat=rbind()
  for(i in inds){
    node=parseGSEAResult$Nodes[[i]]
    es.values=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'ES_PROFILE'),' '))),0)
    hit.index=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'HIT_INDICES'),' '))),nrow(g.rnk))
    m.inds=which.max(abs(es.values))
    es=as.numeric(XML::xmlGetAttr(node,'ES'))
    np=as.numeric(XML::xmlGetAttr(node,'NP'))
    FDR=as.numeric(XML::xmlGetAttr(node,'FDR'))
    nes=as.numeric(XML::xmlGetAttr(node,'NES'))
    title=gsub('^gene_sets.gmt#','',XML::xmlGetAttr(node,'GENESET'))
    length(hit.index)
    all.dat=rbind(all.dat,data.frame(Index=hit.index,ES=es.values,Term=rep(title,length(es.values))
                                     ,Type=c('A',rep('V',length(es.values)-2),'A')))
    all.info=rbind(all.info,c(title,es.values[m.inds[1]],
                              hit.index[m.inds[1]],es,nes,np,FDR))
  }
  all.info=crbind2DataFrame(all.info)
  #all.info
  
  cbPalette=color
  
  all.dat$Colour=cbPalette[as.numeric(as.factor(all.dat$Term))]
  
  col_mp=unique(cbind(as.character(all.dat$Term),all.dat$Colour))[,2]
  names(col_mp)=unique(cbind(as.character(all.dat$Term),all.dat$Colour))[,1]
  glb=unique(unlist(lapply(strsplit(all.info[,1],'_'),function(x){return(x[1])})))
  if(length(glb)==1){
    #g.labels=paste0(gsub(paste0('^',glb,'_'),'',g.labels))
    desc=gsub(paste0('^',glb,'_'),'',all.info[,1])
  }else{
    desc=all.info[,1]
  }
  ndesc=c()
  for(de in desc){
    #de=desc[1]
    if(nchar(de)>50){
      d2=paste0(substr(de,0,47),'...')
      ndesc=c(ndesc,d2)
    }else{
      ndesc=c(ndesc,de)
    }
  }
  g.labels=paste0(ndesc,'\nES=',signif(all.info[,4],2),',NES=',signif(all.info[,5],2),',P=',signif(all.info[,6],2),',FDR=',signif(all.info[,7],2))[match(names(col_mp),all.info[,1])]
  
  #dat=data.frame(Index=hit.index,ES=es.values)
  p=ggplot(data=all.dat, aes(x=Index, y=ES)) +geom_line(aes(colour=Term))+xlim(0,nrow(g.rnk))
  if(length(glb)==1){
    p=p+labs(title=paste0('Enrichment plot ',glb,' terms'))+theme(plot.title = element_text(hjust = 0.5))
  }
  p=p+scale_colour_manual(values=col_mp
                          ,breaks = names(col_mp)
                          ,labels = g.labels
  )
  p=p+ylab('Enrichment score')+theme_bw()
  #p+guides(color = FALSE)
  p=p+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                    , color="grey"
                    ,linetype="dashed")
  
  p=p+theme(
    legend.position='none'
    ,axis.title.x=element_blank()
    ,axis.text.x = element_blank(),axis.ticks.x = element_blank()
    ,plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"))
  #p+geom_label()
  es.min=min(all.dat$ES)
  ymin=es.min-(max(all.dat$ES)-es.min)*0.1
  
  lh=(es.min-ymin)*0.7
  
  dt1=all.dat
  dt2=all.dat
  dt1$Height=rep(ymin,nrow(all.dat))-lh*(as.numeric(as.factor(all.dat$Term))-1)
  dt2$Height=rep(ymin+(es.min-ymin)*0.7,nrow(all.dat))-lh*(as.numeric(as.factor(all.dat$Term))-1)
  dt1=dt1[which(dt1$Type!='A'),]
  dt2=dt2[which(dt2$Type!='A'),]
  #dt=rbind(dt1)
  p1=ggplot()
  #es.text=rbind()
  for(g in unique(dt1$Term)){
    cl=unique(dt1$Colour[which(dt1$Term==g)])[1]
    p1=p1+geom_line(data = rbind(dt1[which(dt1$Term==g),],dt2[which(dt1$Term==g),])
                    ,aes(x=Index,y=Height,group=Index),col=cl)
    #h1=dt1$Height[which(dt1$Term==g)][1]
    #es.text=rbind(es.text,c(g,0,h1))
  }
  #es.text=crbind2DataFrame(es.text)
  #all.info$SX=es.text[match(all.info[,1],es.text[,1]),2]
  #all.info$SY=es.text[match(all.info[,1],es.text[,1]),3]
  
  #p1=p1+geom_text(data = all.info,aes(x=SX,y=SY,label = paste0('ES=',signif(V4,2),',NES=',signif(V5,2),',P=',signif(V6,2),',FDR=',signif(V7,2)))
  #              ,vjust =-0, hjust = 0)
  
  p1=p1+theme_bw()+theme(legend.position='none',axis.title.y=element_blank(),axis.text.y = element_blank()
                         ,axis.title.x=element_blank(),axis.text.x = element_blank()
                         ,axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank()
                         ,axis.line.x.bottom = element_blank()
                         ,plot.margin=unit(c(0, 0.2, 0, 0.1), "inches")
                         ,axis.line = element_blank()
  )
  
  #ggpubr::ggarrange(p,p1, ncol = 1, nrow = 2,heights = c(1,0.1*length(inds)),align = "v")
  
  p2=ggplot(data=data.frame(Risk=c(0,g.rnk$V2,0),Index=c(1,1:nrow(g.rnk),nrow(g.rnk))),aes(y=Risk,x=Index))+geom_line()+theme_bw()
  p2=p2+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                      , color="grey"
                      ,linetype="dashed")
  p2=p2+theme(plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"))+ylab('Rank')+xlab('Rank in ordered dataset')
  p2=p2+geom_text(data=data.frame(Xl=c(0),Yl=c(0)),aes(x=0,y=0,label = left),vjust =1, hjust = 0)+geom_text(data=data.frame(Xl=c(nrow(g.rnk)),Yl=c(0)),
                                                                                                            aes(x=nrow(g.rnk),y=0,label = right),vjust =0, hjust = 1)
  g.h=0.1*length(inds)
  if(g.h>0.8){
    g.h=0.8
  }
  gal=ggpubr::ggarrange(p,p1,p2, ncol = 1, nrow = 3,heights = c(1,g.h,0.6),align = "v",common.legend = TRUE,legend = "right")
  return(gal)
}

plot_GSEA_By_node_wb <- function(parseGSEAResult,
                                 index=1,
                                 col = mycolor,
                                 TermName=NULL,
                                 left=NULL,
                                 right=NULL){
  library(ggplot2)
  if(is.null(parseGSEAResult$TEMPLATE)){
    if(is.null(left)){
      left='RankTop'
    }
    if(is.null(right)){
      right='RankBottom'
    }
  }
  if(is.null(left)){
    left=parseGSEAResult$TEMPLATE[1]
  }
  if(is.null(right)){
    right=parseGSEAResult$TEMPLATE[2]
  }
  if(!is.null(TermName)){
    ind=which(parseGSEAResult$EnrichTable[,1]==TermName)
    if(length(ind)==0){
      print(paste0(TermName,' Not Found!'))
      return(NA)
    }else{
      ind=ind[1]
    }
  }else{
    ind=index
    if(ind>nrow(parseGSEAResult$EnrichTable)){
      print(paste0(ind,' out range!'))
      return(NA)
    }
  }
  node=parseGSEAResult$Nodes[[ind]]
  g.rnk=parseGSEAResult$Rank
  es.values=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'ES_PROFILE'),' '))),0)
  hit.index=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'HIT_INDICES'),' '))),nrow(g.rnk))
  es=as.numeric(XML::xmlGetAttr(node,'ES'))
  nes=as.numeric(XML::xmlGetAttr(node,'NES'))
  np=as.numeric(XML::xmlGetAttr(node,'NP'))
  FDR=as.numeric(XML::xmlGetAttr(node,'FDR'))
  
  dat=data.frame(Index=hit.index,ES=es.values)
  p=ggplot(data=dat, aes(x=Index, y=ES)) +geom_line(aes(colour='darkgreen',size=2))+xlim(0,nrow(g.rnk))+theme_bw()+ylab('Enrichment score')+labs(title=gsub('^gene_sets.gmt#','',XML::xmlGetAttr(node,'GENESET')))+theme(plot.title = element_text(hjust = 0.5))
  p=p+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                    , color="grey"
                    ,linetype="dashed")
  
  #p+geom_text(data=data.frame(Xl=c(0),yl=c(min(es.values))),aes(x=0,y=min(es.values),label = paste0('ES=',signif(es,2),'\nNES=',signif(nes,2)
  #                                                              ,'\nP=',signif(np,2),'\nFDR=',signif(FDR,2)))
  #            ,vjust =0, hjust = 0)
  
  if(es<0){
    p=p+geom_text(data=data.frame(Xl=c(0),yl=c(min(es.values))),aes(x=0,y=min(es.values),label = paste0('ES=',signif(es,2),'\nNES=',signif(nes,2),'\nP=',signif(np,2),'\nFDR=',signif(FDR,2)))
                  ,vjust =0, hjust = 0)
  }else{
    p=p+geom_text(data=data.frame(Xl=c(0),yl=c(min(es.values))),aes(x=nrow(g.rnk),y=max(es.values),label = paste0('ES=',signif(es,2),'\nNES=',signif(nes,2),'\nP=',signif(np,2),'\nFDR=',signif(FDR,2)))
                  ,vjust =1, hjust = 1)
  }
  es.min=min(dat$ES)
  if(es.min>0){
    es.min=0
  }
  
  ymin=es.min-(max(dat$ES)-es.min)*0.1
  dt1=dat
  dt2=dat
  dt1$Height=rep(ymin,nrow(dat))
  dt2$Height=rep(ymin+(es.min-ymin)*0.7,nrow(dat))
  p1=p+geom_line(data = rbind(dt1,dt2),aes(x=Index,y=Height,group=Index))
  p1=p1+ggforce::geom_link(data=data.frame(x=c(0,nrow(g.rnk)),y=c(ymin,ymin),xend=c(nrow(g.rnk)/2,nrow(g.rnk)/2),yend=c(ymin,ymin))
                           ,aes(x=x,y=y,xend=xend,yend=yend
                                ,alpha=1-..index..
                                ,colour=col
                                ,size=50
                           ))
  p1=p1+theme(legend.position='none',axis.title.x=element_blank(),axis.text.x = element_blank())
  p1=p1+geom_text(data=data.frame(Xl=c(0),yl=c(ymin)),
                  aes(x=0,y=ymin,label = left),vjust =0.5, hjust = 0)+
    geom_text(data=data.frame(Xl=c(nrow(g.rnk)),yl=c(ymin)),
              aes(x=nrow(g.rnk),y=ymin,label = right),vjust =0.5, hjust = 1)
  p1=p1+theme(plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"))
  
  p2=ggplot(data=data.frame(Risk=c(0,g.rnk$V2,0),Index=c(1,1:nrow(g.rnk),nrow(g.rnk))),
            aes(y=Risk,x=Index))+geom_line()+theme_bw()
  p2=p2+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                      , color="grey"
                      ,linetype="dashed")
  p2=p2+theme(plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"))+ylab('Rank')+xlab('Rank in ordered dataset')
  gal=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.6),align = "v")
  return(gal)
}

library(ComplexHeatmap)

library(scales)
library(ggsci)
mycolor <- pal_d3(alpha =1)(9)
show_col(mycolor)
dir.create('scripts')

library(data.table)
library(stringr)
dir.create('00_origin_datas')


####
gene_type <- read.delim('GeneTag.genecode.v32.txt', header = T)
gene_type <- gene_type[!duplicated(gene_type$ENSGID), ]
rownames(gene_type) <- gene_type$ENSGID
table(gene_type$TYPE)
gene_type <- crbind2DataFrame(gene_type)
genes_protein <- gene_type[gene_type$TYPE == 'protein_coding', ]$SYMBOL
genes_lncRNA <- gene_type[gene_type$TYPE == 'lncRNA', ]$SYMBOL
str(genes_protein)
str(genes_lncRNA)

# #  ####
# # 
# library(gelnet)
# library(dplyr)
# library(biomaRt)
# install.packages('synapseClient')
# library(synapseClient)
# synapseLogin()
# 
# ### 创建qenes2hugo函数
# #将Ensemble ID转换为HUGO Symbols
# # Use srcType = "ensembl gene id" for Ensembl IDs
# # Use srcType= "entrezgene" for Entrez IDs
# genes2hugo <- function(v,srcType = "ensembl gene id" )
# {
#   ## Retrieve the EMSEMBI -> HUGO mappin
#   gensemb1<- biomaRt::useMart( "ENSEMBI MART ENSEMBL",
#                                host="www.ensembl.org",
#                                dataset="hsapiens gene ensembl")
#   ID <- biomaRt::getBM( attributes=c(srcType,"hgnc symbol"),
#                         filters=srcType,
#                         values=v,mart=ensembl)
#   ## Make sure there was at least one mapping
#   if( nrow(ID) < 1 ) top("No IDs mapped successfully")
#   ## Drop empty duds
#   j <- which( ID[,2] == '')
#   if( length(j) > 0 ) ID <- ID[-j,]
#   stopifnot( all( ID[,1] %in% v ) )
#   ID
# }
# 
# ### main.train 
# ## fnOut - 
# ## fnGenes- 
# main.train <- function( fnOut = "pcbc-stemsiq.tsv",fnGenes = NULL )
# {
#   ## Load RNAseg data
#   synRNA <- synGet("syn2701943",downloadLocation ="/data/PCBC")
#   X <- read.delim( synRNA@filePath ) %>%
#     tibble::column_to_rownames( "tracking_id" ) %>%
#     as.matrix()
#   
#   ## Retrieve metadata
#   synMeta <- synTableQuery( "SELECT UID,Diffname_short FROM syn3156503")
#   Y <- synMeta@values %>%
#     mutate( UID = gsub ("-",".",UID)) %>%
#     tibble::column_to_rownames("UID")
#   
#   ## Retrieve the labels from the metadata
#   y <- Y[colnames (X),]
#   names(y) <- colnames (X)
#   
#   ## Fix the missing labels by hand
#   y["SC11.014BEB.133.5.6.11"] <- "EB"
#   y["SC12.039ECTO.420.436.92.16"] <-"ECTO"
#   
#   ## Drop the splice form ID from the gene names
#   v <- strsplit( rownames (X),"\\.") %>% lapply("[[",1 ) %>% unlist()
#   rownames(X) <- v
#   
#   ## Map Ensembl IDs to HUGO
#   V <- genes2hugo( rownames(X) )
#   X <- X[V[,1],]
#   rownames(X) <- V[,2]
#   
#   ## Reduce the gene set to theprovided list (if applicable)
#   if( is.null( fnGenes ) == FALSE)
#   {
#     vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
#     VE <- genes2hugo( vGenes,"entrezgene")
#     X <- X[intersect( rownames (x), VE[,2] ),]
#   }
#   ## Mean-center the data
#   m <- apply( x,1, mean )
#   X <- X - m
#   
#   ## Identify stem cell samples
#   j <- which( y =="SC")
#   X.tr <- X[,j]
#   X.bk <- X[,-j]
#   
#   ## Train a one-class model
#   mm <- gelnet( t(X.tr),NULL,0,1 )
#   ## Store the signature to a file
#   write.table (mm$w, file = fnOut, sep = "\t",quote = FALSE,col.names = F)
#   
#   ## Perform leave-one-out cross-validation
#   auc <- c()
#   for( i in 1:ncol(X.tr) )
#   {
#     ## Train a model on non-left-out data
#     X1 <- X.tr[,-i]
#     m1 <- gelnet( t(X1), NULL,0,1 )
#     ## Score the left-out sample against the background
#     s.bk <- apply( X.bk, 2, function(z) {cor( ml$w, z, method="sp" )} )
#     s1 <- cor( ml$w,X.tr[,i], method="sp" )
#     ## AUC = P( left-out sample is scored above the backqround )
# 	auc[i] <- sum( sl > s.bk ) / length(s.bk)
# 	cat("Current AUC: ", auc[i],"\n")
# 	cat("Average AUC: ", mean(auc),"\n")
#   }
#   return(auc)
# }
# 
# 
# 
# 
# ###main.predict~
# ## 
# main.predict <- function(fnSig = "pcbc-stemsig.tsv", fnOut = "mRNA_StemScore.tsv")
# {
#   ## Load the signature
#   w <- read.delim( fnSig, header=FALSE,row.names=l ) %>% as.matrix() %>% drop()
# 
#   ## Reduces HUGOIPOSITION gene IDs to just HUGO
#   f <- function(v) unlist( lapply( strsplit( v,"\\|" ),"[[",1) )
#   s <- synGet("syn4976369",downloadLocation ="/data/pancan")
#     X <- read.delim( s@filePath, as.is=TRUE,check.names=FALSE ) %>% ## Read the raw values
#       filter( !grepl( "\\?", gene_id ) ) %>%  ## Drop genes with no mapping to HUGO
#       mutate( gene_id = f( gene_id ) ) %>%  ## Clip gene ids to HUGO
#       filter( gene_id %in% names(w) )  ## Reduce to the signature's gene set
#     
#     ## SIC35E2 has multiple entries with the same HUGO id
#     ## Keep the first entry only
#     j <- grep("SLC35E2",X[,1] )
#     if( length(j) > 1 )
#       X <- X[-j[-1],]
#     
#     ## Convert to a matrix
#     rownames(X) <- NULL
#     X <- X %>% tibble::column_to_rownames( "gene_id") %>% as.matrix()
#     
#     ## Reduce the siqnature to the common set of genes
#     stopifnot( all( rownames(X) %in% names (w) ) )
#     w <- w[rownames (X) ]
#     
#     ####### Score via Spearman correlation
#     s <- apply(X,2,function(z) {cor(z,w,method ="sp",use = "complete.obs")})
#     
#     ## Scale the scores to be between 0 and 1
#     s <- s - min(s)
#     s <- s / max(s)
#     
#     write.table(cbind(s), file = fnOut, sep = "\t", quote = FALSE, col,names = FALSE)
#   }
#   
#   main <- function()
#   {
#     # 
#     main.train()
#     
#     # 
#     main.predict()
#   }














# #####
mRNAsi <- openxlsx::read.xlsx('00_origin_datas/PMID_29625051_mRNAsi.xlsx', sheet = 1,startRow = 1)
tcga_mRNAsi <- mRNAsi[mRNAsi$cancer.type == 'LUAD',]
tcga_mRNAsi <- na.omit(tcga_mRNAsi)
dim(tcga_mRNAsi)  # 569
rownames(tcga_mRNAsi) <- substr(tcga_mRNAsi$TCGAlong.id,1,15)

#  #############
#  #####
tcga_cli2 <- read.delim('00_origin_datas/TCGA/PMC6066282-TCGA-CDR-clinical.txt',header = T, check.names = F)
tcga_cli2 <- subset(tcga_cli2,
                    type == 'LUAD')
table(tcga_cli2$new_tumor_event_type)
tcga_cli2 <- tcga_cli2[, c("bcr_patient_barcode", "type", "OS.time", "OS", 
                           "PFI.time", "PFI", "DFI.time", "DFI", "DSS.time", "DSS")]
colnames(tcga_cli2) <- c("Samples", "type", "OS.time", "OS", 
                         "PFI.time", "PFI", "DFI.time", "DFI", "DSS.time", "DSS")
tcga_cli2$Samples <- paste0(tcga_cli2$Samples, '-01')


tcga_cli <- read.delim('00_origin_datas/TCGA/Clinical BCR XML.merge.txt',
                       header = T, stringsAsFactors = F)
tcga_cli$A0_Samples <- paste0(tcga_cli$A0_Samples, '-01')
tcga_cli <- tcga_cli[, c("A0_Samples", 
                         "age_at_initial_pathologic_diagnosis",
                         "A18_Sex",
                         "A3_T", "A4_N", "A5_M",
                         "A6_Stage")]
colnames(tcga_cli) <- c("Samples",
                        "Age",
                        "Gender",
                        "A3_T", "A4_N", "A5_M",
                        "Stage")

tcga_cli$A3_T <- gsub('[ab]', '', tcga_cli$A3_T)
tcga_cli$A3_T[tcga_cli$A3_T == ''] <- NA
tcga_cli$A3_T[tcga_cli$A3_T == 'TX'] <- NA
tcga_cli$A4_N[tcga_cli$A4_N == ''] <- NA
tcga_cli$A4_N[tcga_cli$A4_N == 'NX'] <- NA
tcga_cli$A5_M <- gsub('[ab]', '', tcga_cli$A5_M)
tcga_cli$A5_M[tcga_cli$A5_M == ''] <- NA
tcga_cli$A5_M[tcga_cli$A5_M == 'MX'] <- NA
tcga_cli$Stage <- gsub('[ABC]', '', tcga_cli$Stage)
tcga_cli$Stage <- gsub('Stage ', '', tcga_cli$Stage)
tcga_cli$Stage[tcga_cli$Stage == ''] <- NA
tcga_cli$Gender <- ifelse(tcga_cli$Gender == 'MALE', 'Male', 'Female')



tcga_cli <- merge(tcga_cli2, tcga_cli, by = 'Samples', all=T)
# tcga_cli <- merge(tcga_cli, tcga_cli1, by = 'Samples', all=T)

rownames(tcga_cli) <- tcga_cli$Samples

tcga_cli <- subset(tcga_cli,
                   !is.na(OS.time) & OS.time > 0)
table(tcga_cli$OS)
# tcga_cli$OS <- ifelse(tcga_cli$OS == 'Alive', 0, 1)
tcga_cli <- crbind2DataFrame(tcga_cli)

colnames(tcga_cli) <- c("Samples","type","OS.time","OS","PFI.time","PFI",
                        "DFI.time","DFI","DSS.time",
                        "DSS","Age","Gender",
                        "T.Stage","N.Stage","M.Stage","Stage")


# 
tcga_exp <- read.table('00_origin_datas/TCGA/TCGA-LUAD-Symbol.txt',
                       sep = '\t', header = T, check.names = F, row.names = 1)
dim(tcga_exp)  # 574

table(substr(colnames(tcga_exp), 14, 17))

tcga_tpm_log2 <- log2(tcga_exp[,] + 1)
tcga_tpm_log2 <- tcga_tpm_log2[rowSums(tcga_tpm_log2) > 0, ]

dim(tcga_tpm_log2)    # 53713 574

tcga_group <- data.frame(row.names = colnames(tcga_tpm_log2),
                         Samples = colnames(tcga_tpm_log2),
                         Groups  = colnames(tcga_tpm_log2),
                         stringsAsFactors = F)
tcga_group$Groups <- substr(colnames(tcga_tpm_log2), 14, 17)
tcga_group <- tcga_group[tcga_group$Groups != '02', ]
table(tcga_group$Groups)
tcga_group$Groups <- ifelse(tcga_group$Groups == '01', 'Tumor', 'Adjacent')

tmr_samples <- intersect(intersect(tcga_cli$Samples, colnames(tcga_tpm_log2)),rownames(tcga_mRNAsi))
length(tmr_samples)  # 494
tcga_cli <- tcga_cli[tmr_samples, ]
dim(tcga_cli)   # 494 16
tcga_mRNAsi <- tcga_mRNAsi[tmr_samples,]

# GSE30219 ###########
GSE30219 <- getGEOExpData('GSE30219')
GSE30219_cli <- GSE30219$Sample
colnames(GSE30219_cli)
GSE30219_cli <- GSE30219_cli[, c("Acc", "Sex", "age", "t-stage", "n-stage",
                                 "status","disease-free survival time")]
colnames(GSE30219_cli) <- c("Samples", "Gender", "Age", "T.Stage", "N.Stage",
                            "OS", "OS.time")

GSE30219_cli$OS<-ifelse(GSE30219_cli$OS == "dead",0 ,1)
GSE30219_cli<-GSE30219_cli[which(GSE30219_cli$OS!= 'NA'),]

GSE30219_cli <- subset(GSE30219_cli,
                       !is.na(OS.time) & OS.time > 0)      
GSE30219_cli$OS.time<-GSE30219_cli$OS.time*365
rownames(GSE30219_cli) <- GSE30219_cli$Samples

GSE30219_anno <- GSE30219$Anno$GPL570
GSE30219_exp <- GSE30219$Exp$GPL570_54675_Data_col1
GSE30219_exp <- exp_probe2symbol_v2(datExpr = GSE30219_exp,
                                    anno = GSE30219_anno[, c(1,11)])
GSE30219_exp <- GSE30219_exp[str_split_fixed(rownames(GSE30219_exp), ' /// ', 2)[, 2] == '', ]
GSE30219_com_samples <- intersect(GSE30219_cli$Samples, colnames(GSE30219_exp))
GSE30219_cli <- GSE30219_cli[GSE30219_com_samples, ]

range(GSE30219_exp)
dim(GSE30219_exp)

# GSE19188 ###########
GSE19188 <- getGEOExpData('GSE19188')
GSE19188_cli <- GSE19188$Sample
colnames(GSE19188_cli)
GSE19188_cli <- GSE19188_cli[, c("Acc", "status","overall survival")]
colnames(GSE19188_cli) <- c("Samples", "OS", "OS.time")

GSE19188_cli<-GSE19188_cli[which(GSE19188_cli$OS!= 'Not available'),]
GSE19188_cli$OS<-ifelse(GSE19188_cli$OS == "alive",1 ,0)


GSE19188_cli <- subset(GSE19188_cli,
                       !is.na(OS.time) & OS.time > 0)   
GSE19188_cli$OS.time<-as.numeric(GSE19188_cli$OS.time)
GSE19188_cli$OS.time<-GSE19188_cli$OS.time*30
rownames(GSE19188_cli) <- GSE19188_cli$Samples

GSE19188_anno <- GSE19188$Anno$GPL570
GSE19188_exp <- GSE19188$Exp$GPL570_54675_Data_col1
GSE19188_exp <- exp_probe2symbol_v2(datExpr = GSE19188_exp,
                                    anno = GSE19188_anno[, c(1,11)])

GSE19188_exp <- GSE19188_exp[str_split_fixed(rownames(GSE19188_exp), ' /// ', 2)[, 2] == '', ]
GSE19188_com_samples <- intersect(GSE19188_cli$Samples, colnames(GSE19188_exp))
GSE19188_cli <- GSE19188_cli[GSE19188_com_samples, ]

range(GSE19188_exp)
dim(GSE19188_exp)





#  ####
dir.create('01_miRNAsi')

tcga_mRNAsi_cli <- tcga_cli
tcga_mRNAsi_cli$mRNAsi <- tcga_mRNAsi$mRNAsi
tcga_mRNAsi_cli$Age <- ifelse(tcga_mRNAsi_cli$Age >60,'>60','<=60')
colnames(tcga_mRNAsi_cli)

mRNAsi_Stage <- data.frame(tcga_mRNAsi_cli[, c('Stage', "mRNAsi")])
mRNAsi_Stage <- na.omit(mRNAsi_Stage)
mRNAsi_Stage_violin <- mg_violin(mRNAsi_Stage, melt=TRUE, ylab='mRNAsi',
                                  leg.title='', 
                                  xlab = 'Stage',
                                  jitter = F,
                                  legend.pos='bl',
                                  show_compare = F)

mRNAsi_Stage_violin   # 0.064

mRNAsi_T.Stage <- data.frame(tcga_mRNAsi_cli[, c('T.Stage', "mRNAsi")])
mRNAsi_T.Stage <- na.omit(mRNAsi_T.Stage)
mRNAsi_T.Stage_violin <- mg_violin(mRNAsi_T.Stage, melt=TRUE, ylab='mRNAsi',
                                 leg.title='', 
                                 xlab = 'T.Stage',
                                 jitter = F,
                                 legend.pos='bl',
                                 show_compare = F)

mRNAsi_T.Stage_violin   # 0.019

mRNAsi_N.Stage <- data.frame(tcga_mRNAsi_cli[, c('N.Stage', "mRNAsi")])
mRNAsi_N.Stage <- na.omit(mRNAsi_N.Stage)
mRNAsi_N.Stage_violin <- mg_violin(mRNAsi_N.Stage, melt=TRUE, ylab='mRNAsi',
                                 leg.title='', 
                                 xlab = 'N.Stage',
                                 jitter = F,
                                 legend.pos='bl',
                                 show_compare = F)

mRNAsi_N.Stage_violin   # 0.29

mRNAsi_M.Stage <- data.frame(tcga_mRNAsi_cli[, c('M.Stage', "mRNAsi")])
mRNAsi_M.Stage <- na.omit(mRNAsi_M.Stage)
mRNAsi_M.Stage_violin <- mg_violin(mRNAsi_M.Stage, melt=TRUE, ylab='mRNAsi',
                                 leg.title='', 
                                 xlab = 'M.Stage',
                                 jitter = F,
                                 legend.pos='bl',
                                 show_compare = F)

mRNAsi_M.Stage_violin   # 0.076

mRNAsi_Gender <- data.frame(tcga_mRNAsi_cli[, c('Gender', "mRNAsi")])
mRNAsi_Gender <- na.omit(mRNAsi_Gender)
mRNAsi_Gender_violin <- mg_violin(mRNAsi_Gender, melt=TRUE, ylab='mRNAsi',
                               leg.title='', 
                               xlab = 'Gender',
                               jitter = F,
                               legend.pos='bl',
                               show_compare = F)

mRNAsi_Gender_violin  # 0.00027


library(survcomp)
library(survminer)
tcga_model_data <- cbind(tcga_mRNAsi_cli[, c("OS.time", "OS")],
                         t(tcga_tpm_log2[,tcga_mRNAsi_cli$Samples]))
tra.data$mRNAsi <- tcga_mRNAsi$mRNAsi
tra.data.point <- surv_cutpoint(tra.data, time = "OS.time", event = "OS",
                                variables = 'mRNAsi')
tra.cutoff <- as.numeric(summary(tra.data.point)[1])
tra.cutoff

tr.km <- ggplotKMCox(data.frame(tra.data$OS.time / 365,
                                tra.data$OS,
                                ifelse(tcga_mRNAsi$mRNAsi>=tra.cutoff,'High','Low')),
                     title = 'mRNAsi',
                     labs = c('High','Low'))
tr.km

fig1 <- cowplot::plot_grid(mRNAsi_Gender_violin,
                           mRNAsi_Stage_violin,
                           mRNAsi_T.Stage_violin,
                           mRNAsi_N.Stage_violin,
                           mRNAsi_M.Stage_violin,
                           tr.km,
                   ncol = 3,nrow = 2,
                   labels=LETTERS[1:6],
                   label_fontfamily = 'Times',
                   label_size = 14)
fig1
ggsave(plot = fig1,
       filename = '01_miRNAsi/fig1.pdf',
       width = 12, height = 10)
#  #####

tcga_limma_dat <- tcga_tpm_log2[rownames(tcga_tpm_log2) %in% genes_protein, ]     #
tcga_limma_dat <- tcga_limma_dat[rowSums(tcga_limma_dat) > 0, ]

tcga_mRNAsi_cli$Group <- ifelse(tcga_mRNAsi$mRNAsi>=tra.cutoff,'High','Low')
table(tcga_mRNAsi_cli$Group)

tcga_mRNAsi_limma <- mg_limma_DEG(exp = tcga_limma_dat[, tcga_mRNAsi_cli$Samples],
                                group = tcga_mRNAsi_cli$Group,
                                ulab = 'High',
                                dlab = 'Low')
tcga_mRNAsi_limma_res <- tcga_mRNAsi_limma$DEG
tcga_mRNAsi_limma_res <- na.omit(tcga_mRNAsi_limma_res)
head(tcga_mRNAsi_limma_res)
tcga_mRNAsi_limma_filtered <- subset(tcga_mRNAsi_limma_res,
                                   abs(logFC) > log2(1.5) & adj.P.Val < 0.05)          #
dim(tcga_mRNAsi_limma_filtered)   # 1929
table(tcga_mRNAsi_limma_filtered$logFC > 0)  # 527 1402

tcga_mRNAsi_limma_genes <- rownames(tcga_mRNAsi_limma_filtered)

tcga_mRNAsi_volcano <- mg_volcano_wb(logfc = tcga_mRNAsi_limma_res$logFC,
                                   pvalue = tcga_mRNAsi_limma_res$adj.P.Val,
                                   cutFC = log2(1.5),
                                   cutPvalue = 0.05,
                                   colors = c(mycolor[1], 'grey', mycolor[2]),
                                   legend.pos = 'tl',
                                   leg = 'TCGA')
tcga_mRNAsi_volcano

ggsave(plot = tcga_mRNAsi_volcano,
       filename = '02_WGCNA/fig2A.pdf',
       width = 4, height = 5)

#
tcga_mRNAsi_go_kegg <- enrichmentORA(tcga_mRNAsi_limma_genes,
                                   mp_dbs=c('pathway_KEGG',
                                            'geneontology_Biological_Process',
                                            'geneontology_Cellular_Component',
                                            'geneontology_Molecular_Function'))
tcga_mRNAsi_go_kegg_filtered <- tcga_mRNAsi_go_kegg[tcga_mRNAsi_go_kegg$FDR < 0.05, ]
table(tcga_mRNAsi_go_kegg_filtered$DB)

pdf('02_WGCNA/tcga_mRNAsi_go_kegg.pdf', width = 11, height = 8)
dotplot_batch(tcga_mRNAsi_go_kegg_filtered, 
              dbs =c('geneontology_Biological_Process',
                     'geneontology_Cellular_Component',
                     'geneontology_Molecular_Function',
                     'pathway_KEGG'),top=10,FDR = T)
dev.off()

# WGCNA ####################
dir.create('02_WGCNA')
diff.genes <- tcga_mRNAsi_limma_genes
tcga_tmr_tpm <- tcga_tpm_log2[diff.genes,tcga_mRNAsi_cli$Samples]
tcga_tmr_filter <- tcga_tmr_tpm[which(apply(tcga_tmr_tpm,1,function(x){return(sum(x>=1))})>0.5*ncol(tcga_tpm_log2)),]
library(WGCNA)
tcga_tpm_wgcna <- tcga_tmr_filter
sds <- apply(tcga_tpm_wgcna,1,sd)
tcga_tpm_wgcna <- tcga_tpm_wgcna[which(sds>0),]
dim(tcga_tpm_wgcna)   # 1784 494

tcga_edata.wgcna=as.data.frame(t(tcga_tpm_wgcna))

pdf('02_WGCNA/fig2ABC.pdf',width = 8,height = 8)
power=mg_wgcna_get_power(tcga_edata.wgcna)
dev.off()

power$cutPower   # 5
net=mg_WGCNA_getModule(tcga_edata.wgcna,power = power$cutPower
                       , deepSplit = 2, mergeCutHeight = 0.25
                       , minModuleSize = 30)
length(table(net$Modules[,2]))  # 6


pdf('02_WGCNA/fig2D.pdf',height = 4,width = 6)
plotDendroAndColors(net$Tree, net$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

fig2E=mg_barplot_point(labels = names(table(net$Modules[,2]))
                       ,values = as.numeric(table(net$Modules[,2]))
                       ,point_sizes = 2
                       ,point_cols = names(table(net$Modules[,2]))
                       ,xlab = 'Number of Genes',legend.pos = NULL)
fig2E

savePDF('02_WGCNA/Fig2E.pdf',fig2E,height = 6,width = 6)

#### 
# Calculate eigengenes
MEs = net$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result

# pdf('02_WGCNA/Fig2F.pdf',height = 6,width = 10,onefile = T)
# plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
# dev.off()

#  ####
### 
infil.t <- tcga_mRNAsi$mRNAsi
datTraits <- infil.t
####### Calculate module eigengenes
MEs<-net$MEs
dim(MEs)
## Define numbers of genes and samples
nGenes = ncol(tcga_edata.wgcna)     # 6272
nSamples = nrow(tcga_edata.wgcna)   # 500
## 
modTraitCor = WGCNA::cor(MEs[,rownames(MEDiss)[METree$order]]
                         , datTraits
                         , use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, nSamples)
## 
textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

dev.off()
colnames(modTraitCor) <- 'mRNAsi'
pdf('02_WGCNA/Fig2F.pdf',width = 6,height = 6)
labeledHeatmap(Matrix = data.frame(modTraitCor), 
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = rownames(modTraitCor), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix), setStdMargins = T,
               cex.text = 1, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#
geneModuleMembership <- as.data.frame(signedKME(tcga_edata.wgcna, data.frame(net$MEs), outputColumnName = ""))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

#
all(rownames(tcga_edata.wgcna)==rownames(datTraits))
geneTraitSignificance <- as.data.frame(cor(tcga_edata.wgcna, datTraits, use = 'pairwise.complete.obs'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

###
### Intramodular analysis: identifying genes with high GS and MM
############
get_module_hub_genes=function(net=NULL,module=NULL,trait=NULL,output=NULL,MM=0.6,GS=0.5,pval=0.05){
  modNames = substring(names(net$MEs), 3)
  moduleColors = unname(net$Modules[,2])
  #########
  module = module
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  pdf(output,height = 6,width = 6)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, trait]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",trait),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                     , col = module,lwd=2)
  abline(v = MM, col = "red", lwd = 2, lty = 2)
  abline(h = GS, col = "red", lwd = 2, lty = 2)
  dev.off()
  print(colnames(geneModuleMembership)[column])
  inds=(abs(geneModuleMembership[moduleGenes, column])>MM) & (abs(geneTraitSignificance[moduleGenes, trait])>GS)
  hub.genes=rownames(geneModuleMembership)[moduleGenes][inds]
  return(hub.genes)
}

table(net$Modules[,2])
colnames(geneTraitSignificance)
# module-turquoise  hub基因
turquoise.hub.genes=get_module_hub_genes(net = net,module = "turquoise",trait = "V1",output=NULL
                                     ,MM=0.4,GS=0.4,pval=0.05)
length(turquoise.hub.genes)   # 451
write.table(turquoise.hub.genes,file = 'results/turquoise_hub_genes.txt',sep = "\t",quote = F,row.names = F)

# module-turquoise  模块基因
turquoise.genes=get_module_hub_genes(net = net,module = "turquoise",trait = "V1"
                                 ,MM=0,GS=0,pval=1)
length(turquoise.genes)   # 775

# 
modNames = substring(names(MEs), 3)
moduleColors = unname(net$Modules[,2])
table(moduleColors)
module = "turquoise"
column = match(module, modNames)
moduleGenes = moduleColors==module
table(moduleGenes)
trait="V1"
table(moduleColors)[module]
pdf('02_WGCNA/fig2G.pdf',width = 6,height = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, trait]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ",trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
abline(v = 0.4, col = mycolor[2], lwd = 2, lty = 2)
abline(h = 0.4, col = mycolor[2], lwd = 2, lty = 2)
dev.off()
#  ####################
dir.create('03_model')

tcga_model_data <- cbind(tcga_mRNAsi_cli[, c("OS.time", "OS")],
                         t(tcga_tpm_log2[turquoise.hub.genes, tcga_mRNAsi_cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
tcga_model_data <- na.omit(tcga_model_data)
tcga_model_data <- crbind2DataFrame(tcga_model_data)


coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
#
tra.cox <- t(apply(tcga_model_data[,3:c(ncol(tcga_model_data))],2,function(x){
  vl=as.numeric(x)
  tm=tcga_model_data$OS.time
  ev=tcga_model_data$OS
  #ev=ifelse(ev=='Alive',0,1)
  dat=data.frame(tm,ev,vl)[which(tm > 0 & !is.na(vl)),]
  return(coxFun(dat))
}))
colnames(tra.cox)=c('p.value','HR','Low 95%CI','High 95%CI')
length(which(tra.cox[,1]<0.05))
tra.cox <- na.omit(tra.cox)

dim(tra.cox)

tra.cox1 <- as.data.frame(tra.cox)
# tra.cox1 <- tra.cox1[tra.cox1$HR < 10, ]
# tra.cox1 <- tra.cox1[tra.cox1$HR > 1e-5, ]
dim(tra.cox1)
# pdf('03_model/sig.cox.pdf', width = 5, height = 5)
mg_volcano_risk_wb <- function(logfc,
                               pvalue,
                               symbol=NULL,
                               cutFC=1,
                               cutPvalue=0.05
                               ,showText=NULL
                               ,colors=c(mg_colors[2],'grey',mg_colors[1])
                               ,xlim=NULL,ylim=NULL
                               ,legend.pos='right'
                               ,ylab='-log10(FDR)',
                               leg='State',
                               xlab='log2(FoldChange)'){
  library(ggplot2)
  pos=c(0,0)
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else{
    pos='right'
  }
  cange=rep('None',length(logfc))
  cange[which(logfc>cutFC&pvalue<cutPvalue)]='Risk'
  cange[which(logfc< -cutFC&pvalue<cutPvalue)]='Protect'
  if(is.null(symbol)){
    symbol=rep('',length(logfc))
    showText=NULL
  }
  vo.input=data.frame(logFC=logfc,FDR=pvalue,change=cange,SYMBOL=symbol)
  #print(head(vo.input))
  p1 <- ggplot(data = vo.input, 
               aes(x = logFC, 
                   y = -log10(FDR)))
  p1=p1+geom_point(alpha=0.75, size=2.5, aes(color=change))  
  p1=p1+scale_color_manual(values=colors,limits = c("Protect",'None', "Risk"),name=leg) 
  p1=p1+geom_vline(xintercept=c(-cutFC,cutFC),lty=4,col="black",lwd=0.8)  
  p1=p1+geom_hline(yintercept = -log10(cutPvalue),lty=4,col="black",lwd=0.8)  
  p1=p1+ylab(ylab)+xlab(xlab)
  p1=p1+theme_bw()
  p1=p1+theme(
    axis.text.y=element_text(family="Times",face="plain"), #
    axis.title.y=element_text(family="Times",face="plain"), #
    legend.text=element_text(face="plain", family="Times", colour="black"  #
    ),
    legend.title=element_text(face="plain", family="Times", colour="black" #
    ),
    legend.justification=pos, legend.position=pos
    ,legend.background = element_rect(fill = NA, colour = NA)
  )
  if(is.null(showText)|is.null(symbol)){
    showText=c()
  }
  
  if(length(showText)>0){
    for_label <-vo.input[match(intersect(showText,vo.input$SYMBOL),vo.input$SYMBOL),]
    p1=p1+geom_point(size = 3, shape = 1, data = for_label)+ggrepel::geom_label_repel(
      aes(label = SYMBOL),
      data = for_label,
      color="black"
    )
  }
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(!is.null(xlim)){
    p1=p1+xlim(xlim)
  }
  
  return(p1)
}

filter_genes <- rownames(tra.cox[tra.cox[,1]<0.05, ])
length(filter_genes)  # 345

tra.cox.filtered <- as.data.frame(tra.cox[tra.cox[,1]<0.05, ])
table(tra.cox.filtered$HR > 1)

library(glmnet)
fit1=glmnet(as.matrix(tcga_model_data[,filter_genes])
            #,factor(samps)
            ,cbind(time=tcga_model_data$OS.time,
                   status=tcga_model_data$OS)
            ,family="cox"
            #,family="binomial"
            #,type.measure="deviance"
            ,nlambda=100
            , alpha=1)
set.seed(2023)
cv.fit<-cv.glmnet(as.matrix(tcga_model_data[,filter_genes])
                  #,factor(samps)
                  ,cbind(time=tcga_model_data$OS.time,
                         status=tcga_model_data$OS)
                  ,family="cox"
                  #,family="binomial"
                  #,type.measure="deviance"
                  ,nlambda=100
                  , alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
length(sig.coef)  # 8
cv.fit$lambda.min

lasso.pdf <- mg_plot_lasso(fit1,
                           cv.fit,
                           # lambda = cv.fit$lambda.min,
                           show_text=T,
                           figLabels=c('A','B'))
lasso.pdf
ggsave(plot = lasso.pdf,
       filename = '03_model/fig3AB.pdf',
       width = 8, height = 4)

lasso_genes <- names(sig.coef)

# 
tcga_dat1 <- cbind(time=tcga_model_data$OS.time,
                   status=tcga_model_data$OS,
                   tcga_model_data[, lasso_genes])

fmla <- as.formula(paste0("Surv(time, status) ~"
                          ,paste0(lasso_genes, collapse = '+')))


cox <- coxph(fmla, data =as.data.frame(tcga_dat1))

#cox1 <- step(cox, trace = 0)
lan <- coef(cox)
round(lan, 3)
genes <- names(cox$coefficients)
tra.cox[genes,]
summary(cox1)
# pdf('03_model/fig3C.pdf', width = 6, height = 4, onefile = F)
# survminer::ggforest(cox,data=tcga_dat1)
# dev.off()


lst.modl=createCoxModel(as.matrix(tcga_dat1[, lasso_genes])
                        ,time = tcga_dat1$time
                        ,event = tcga_dat1$status,
                        isStep = F)
lst.modl$Genes
lst.modl$fmla

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)
library(dplyr)
library(ggsci)
fig3C=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  scale_fill_d3() +
  # coord_flip() +
  labs(x = "") +
  labs(y = "LASSO Cox coefficient") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position = 'top')
fig3C
ggsave(plot = fig3C,
       filename = '03_model/fig3C.pdf',
       width = 5, height = 4)


fig3ABC <- cowplot::plot_grid(lasso.pdf,
                           fig3C,
                           labels=c('','C'),
                           label_fontfamily = 'Times',
                           label_size = 14,
                           rel_widths = c(1,0.6))
fig3ABC
ggsave(plot = fig3ABC,
       filename = '03_model/fig3ABC.pdf',
       width = 14, height = 4)

library(survcomp)
library(survminer)
# 
risk.tr <- as.numeric(lan%*%as.matrix(t(tra.data[,genes])))
tcga_model_data$RS <- risk.tr
# tcga_model_data.point <- surv_cutpoint(tcga_model_data, time = "OS.time", event = "OS",
#                                 variables = 'RS')
# tcga.cutoff <- as.numeric(summary(tcga_model_data.point)[1])
tcga.cutoff <- median(tcga_model_data$RS)
tcga_model_data$RiskType <- ifelse(tcga_model_data$RS>=tcga.cutoff,'High','Low')


tcga.roc <- ggplotTimeROC(tcga_model_data$OS.time / 365,
                        tcga_model_data$OS,
                        risk.tr,
                        mks = c(1,3,5))
tcga.km <- ggplotKMCox(data.frame(tcga_model_data$OS.time / 365,
                                tcga_model_data$OS,
                                tcga_model_data$RiskType),
                     title = 'TCGA OS',
                     labs = c('High','Low'))

tcga.roc.km <- cowplot::plot_grid(tcga.km,
                                tcga.roc,
                                nrow = 2)
tcga.roc.km


plotRiskScoreModel=function(riskScore,dat,time,event,cutoff,hetTitle='z-score of expression',hetColor=c('green','black','red')){
  srt.inds=order(riskScore)
  dat=dat[srt.inds,]
  time=time[srt.inds]
  event=event[srt.inds]
  riskScore=riskScore[srt.inds]
  library(ggplot2)
  dt1=data.frame(V1=1:length(riskScore),V2=riskScore,RiskType=ifelse(riskScore>cutoff,'High','Low')) 
  p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +geom_point(stat = 'identity', position = 'dodge')+ggsci::scale_fill_d3()+theme_bw()
  p1=p1+ylab('RiskScore')+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
                                ,axis.title.x=element_blank(),legend.position=c(1,0), legend.justification=c(1,0)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,plot.margin=unit(c(0.1, 0.1, 0, 0.1), "inches")
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  
  dt2=data.frame(V1=1:length(riskScore),V2=time,Status=ifelse(event==1,'Dead','Alive'))  
  p2=ggplot(dt2, aes(x = V1, y = V2, colour = Status,shape =Status)) +geom_point()+ggsci::scale_fill_d3()+theme_bw()
  p2=p2+ylab('Time')+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
                           ,axis.title.x=element_blank(),legend.position=c(1,1), legend.justification=c(1,1)
                           ,legend.background = element_rect(fill = NA, colour = NA)
                           ,plot.margin=unit(c(0, 0.1, 0, 0.1), "inches")
                           ,legend.title = element_text(family="Times",face="plain")
                           ,legend.text = element_text(family="Times",face="plain"))
  
  data=as.data.frame(scale(dat))
  hc.r = hclust(dist(t(data)))
  data=data[,hc.r$order]
  data$ID <- 1:nrow(dat)
  #colnames(data)
  data_m <- reshape2::melt(data, id.vars=c("ID"))
  colnames(data_m)=c('ID','V1','V2')
  data_m$V2[which(data_m$V2>mean(data_m$V2)+3*sd(data_m$V2))]=mean(data_m$V2)+3*sd(data_m$V2)
  data_m$V2[which(data_m$V2<mean(data_m$V2)-3*sd(data_m$V2))]=mean(data_m$V2)-3*sd(data_m$V2)
  
  data_m$V1=mg_str_outline(data_m$V1,isCut = T,n=50)
  #print(data_m$V1)
  #print(head(data_m))
  #data_m[1:20,]
  p3 <- ggplot(data_m, aes(x=ID,y=V1)) 
  p3 <- p3 + geom_tile(aes(fill=V2))
  p3=p3+scale_fill_gradient2(low = hetColor[1],mid=hetColor[2], high = hetColor[3])
  p3=p3+theme_bw()
  p3=p3+ labs(fill=hetTitle) 
  #p3=p3+guides(fill=guide_legend(title="New Legend Title"))
  p3=p3+xlab('Samples')
  #p3 <- p3 + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  p3=p3+theme(axis.text.y=element_text(family="Times",face="plain")
              ,axis.text.x=element_blank()
              #,axis.title.x=element_blank()
              ,axis.title.y=element_blank()
              ,legend.position='bottom'
              #,legend.justification=c(1,1)
              #,legend.background = element_rect(fill = NA, colour = NA)
              ,plot.margin=unit(c(0, 0.1, 0.1, 0.1), "inches")
              ,legend.title = element_text(family="Times",face="plain")
              ,legend.text = element_text(family="Times",face="plain"))
  
  g1=ggpubr::ggarrange(p1,p2,p3, ncol = 1, nrow = 3,heights = c(0.5,0.5,1),align = "v")
  return(g1)
}

tcga_genes_data<-tcga_model_data[,genes]
tcga_genes_data$OS <- tcga_model_data$OS
tcga_genes_data$OS.time <- tcga_model_data$OS.time

dim(tcga_genes_data)
tcga_point <- plotRiskScoreModel(riskScore = risk.tr,
                                       dat = tcga_genes_data[,genes],
                                       time = tcga_genes_data$OS.time,
                                       event = tcga_genes_data$OS,
                                       cutoff = tcga.cutoff)

tcga_point
fig3DE <- cowplot::plot_grid(tcga.roc.km,
                           tcga_point,
                           ncol = 2,
                           labels=c('D','E'),
                           label_fontfamily = 'Times',
                           label_size = 14)
fig3DE
ggsave(plot = fig3DE,
       filename = '03_model/fig3DE.pdf',
       width = 14, height = 8)

tcga_RS_mRNAsi_cor <- cor_point(x = tcga_mRNAsi$mRNAsi,
                             y = tcga_model_data$RS,
                             top_col = mycolor[1],
                             right_col = mycolor[2],
                             xlab = 'mRNAsi',
                             ylab = 'RiskScore')
tcga_RS_mRNAsi_cor
ggsave(plot = tcga_RS_mRNAsi_cor,
       filename = '03_model/tcga_RS_mRNAsi_cor.pdf',
       width = 6, height = 4)

#  GSE19188
GSE19188_model_data <- cbind(GSE19188_cli[, c("OS.time", "OS")],
                             t(GSE19188_exp[turquoise.hub.genes, GSE19188_cli$Samples]))
colnames(GSE19188_model_data) <- gsub('-', '_', colnames(GSE19188_model_data))
GSE19188_model_data <- crbind2DataFrame(GSE19188_model_data)
GSE19188.genes <- intersect(genes, colnames(GSE19188_model_data))
fmla.GSE19188 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                   ,paste0(GSE19188.genes,collapse = '+')))
cox.GSE19188 <- coxph(fmla.GSE19188, data =as.data.frame(GSE19188_model_data))
GSE19188_lan <- coef(cox.GSE19188)

setdiff(names(lan),colnames(GSE19188_model_data))
risk.GSE19188=as.numeric(GSE19188_lan%*%as.matrix((t(GSE19188_model_data[,names(GSE19188_lan)]))))
risk.GSE19188

library(survminer)
# risk.GSE19188=as.numeric(lan%*%as.matrix(t(GSE19188_model_data[,names(lan)])))

GSE19188_model_data$RS <- risk.GSE19188
# GSE19188_model_data$RS <- predict(cox.GSE19188)
GSE19188.data.point <- surv_cutpoint(GSE19188_model_data, time = "OS.time", event = "OS",
                                     variables = 'RS')
GSE19188.cutoff <- as.numeric(summary(GSE19188.data.point)[1])
GSE19188.cutoff <- median(risk.GSE19188)



GSE19188.os.roc <- ggplotTimeROC(GSE19188_model_data$OS.time/365,
                                 GSE19188_model_data$OS,
                                 risk.GSE19188,
                                 mks = c(1,3,5))
GSE19188.os.roc
library(survival)
GSE19188.os.km <- ggplotKMCox(data.frame(GSE19188_model_data$OS.time/365 ,
                                         GSE19188_model_data$OS,
                                         ifelse(risk.GSE19188>=GSE19188.cutoff,'High','Low')),
                              title = 'GSE19188',
                              labs = c('High','Low'))
GSE19188.os.km

GSE19188.roc.km <- cowplot::plot_grid(GSE19188.os.km,
                                      GSE19188.os.roc,
                                      ncol = 1,
                                      labels=c('F'))
GSE19188.roc.km

#  GSE30219
GSE30219_model_data <- cbind(GSE30219_cli[, c("OS.time", "OS")],
                             t(GSE30219_exp[turquoise.hub.genes, GSE30219_cli$Samples]))
colnames(GSE30219_model_data) <- gsub('-', '_', colnames(GSE30219_model_data))
GSE30219_model_data <- crbind2DataFrame(GSE30219_model_data)
GSE30219.genes <- intersect(genes, colnames(GSE30219_model_data))
fmla.GSE30219 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                   ,paste0(GSE30219.genes,collapse = '+')))
cox.GSE30219 <- coxph(fmla.GSE30219, data =as.data.frame(GSE30219_model_data))
GSE30219_lan <- coef(cox.GSE30219)

setdiff(names(lan),colnames(GSE30219_model_data))
risk.GSE30219=as.numeric(GSE30219_lan%*%as.matrix((t(GSE30219_model_data[,names(GSE30219_lan)]))))
risk.GSE30219 

library(survminer)
# risk.GSE30219=as.numeric(lan%*%as.matrix(t(GSE30219_model_data[,names(lan)])))

GSE30219_model_data$RS <- risk.GSE30219
# GSE30219_model_data$RS <- predict(cox.GSE30219)
GSE30219.data.point <- surv_cutpoint(GSE30219_model_data, time = "OS.time", event = "OS",
                                     variables = 'RS')
GSE30219.cutoff <- as.numeric(summary(GSE30219.data.point)[1])
GSE30219.cutoff <- median(risk.GSE30219)



GSE30219.os.roc <- ggplotTimeROC(GSE30219_model_data$OS.time/365,
                                 GSE30219_model_data$OS,
                                 risk.GSE30219,
                                 mks = c(1,3,5))
GSE30219.os.roc
library(survival)
GSE30219.os.km <- ggplotKMCox(data.frame(GSE30219_model_data$OS.time/365 ,
                                         GSE30219_model_data$OS,
                                         ifelse(risk.GSE30219>=GSE30219.cutoff,'High','Low')),
                              title = 'GSE30219',
                              labs = c('High','Low'))
GSE30219.os.km

GSE30219.roc.km <- cowplot::plot_grid(GSE30219.os.km,
                                      GSE30219.os.roc,
                                      ncol = 1,
                                      labels=c('G'))
GSE30219.roc.km
fig3FG <- cowplot::plot_grid(GSE19188.roc.km,
                                      GSE30219.roc.km,
                                      ncol = 2,
                                      label_fontfamily = 'Times',
                                      label_size = 14)

fig3FG
ggsave(plot = fig3FG,
       filename = '03_model/fig3FG.pdf',
       width = 14, height = 8)

#  #####
dir.create('04_RS_diff.genes')

tcga_limma_dat <- tcga_tpm_log2[rownames(tcga_tpm_log2) %in% genes_protein, ]     #
tcga_limma_dat <- tcga_limma_dat[rowSums(tcga_limma_dat) > 0, ]

table(tcga_model_data$RiskType)

tcga_risk_limma <- mg_limma_DEG(exp = tcga_limma_dat[, rownames(tcga_model_data)],
                                group = tcga_model_data$RiskType,
                                ulab = 'High',
                                dlab = 'Low')
tcga_risk_limma_res <- tcga_risk_limma$DEG
tcga_risk_limma_res <- na.omit(tcga_risk_limma_res)
head(tcga_risk_limma_res)
tcga_risk_limma_filtered <- subset(tcga_risk_limma_res,
                                   abs(logFC) > 1 & adj.P.Val < 0.05)          #
dim(tcga_risk_limma_filtered)    # 400
table(tcga_risk_limma_filtered$logFC > 0)   # 221 179

write.table(tcga_risk_limma_filtered,file = 'results/tcga_risk_limma_filtered.txt',sep = "\t",quote = F,row.names = F)

tcga_risk_limma_genes <- rownames(tcga_risk_limma_filtered)

tcga_risk_volcano <- mg_volcano_wb(logfc = tcga_risk_limma_res$logFC,
                                   pvalue = tcga_risk_limma_res$adj.P.Val,
                                   cutFC = 1.5,
                                   cutPvalue = 0.05,
                                   colors = c(mycolor[1], 'grey', mycolor[2]),
                                   legend.pos = 'tl',
                                   leg = 'TCGA')
tcga_risk_volcano
ggsave(plot = tcga_risk_volcano,
       filename = '04_RS_diff.genes/fig4A.pdf',
       width = 5, height = 6)
#
tcga_risk_go_kegg <- enrichmentORA(tcga_risk_limma_genes,
                                   mp_dbs=c('pathway_KEGG',
                                            'geneontology_Biological_Process',
                                            'geneontology_Cellular_Component',
                                            'geneontology_Molecular_Function'))
tcga_risk_go_kegg_filtered <- tcga_risk_go_kegg[tcga_risk_go_kegg$FDR < 0.05, ]
table(tcga_risk_go_kegg_filtered$DB)

pdf('04_RS_diff.genes/fig4BC.pdf', width = 10, height = 6)
dotplot_batch(tcga_risk_go_kegg_filtered, 
              dbs =c('geneontology_Biological_Process',
                     'geneontology_Cellular_Component'
                     #'geneontology_Molecular_Function',
                     #'pathway_KEGG'
                     ),top=10,FDR = T)
dev.off()


#  ###############
dir.create('05_RS_immune')

tcga_RS_cli <- tcga_cli
tcga_RS_cli$RiskScore <- tcga_model_data$RS
tcga_RS_cli$RiskType <- tcga_model_data$RiskType

tme1 <- openxlsx::read.xlsx('00_origin_datas/PMID.28052254.xlsx', sheet = 1)
tme1.list1 <- split(tme1$Metagene, tme1$Cell.type)
tme1.list2 <- split(tme1$Metagene, tme1$Immunity)

tme2 <-clusterProfiler::read.gmt('00_origin_datas/29_immu_signature_pmid_30594216.txt')
tme2 <- tme2[tme2$gene != '', ]
tme2.list <- split(tme2$gene, tme2$ont)

ImmuCellAI <- clusterProfiler::read.gmt('00_origin_datas/ImmuCellAI_marker_genes.txt')
ImmuCellAI <- ImmuCellAI[ImmuCellAI$gene != '', ]
table(ImmuCellAI$ont)
ImmuCellAI.list <- split(ImmuCellAI$gene, ImmuCellAI$ont)



save(tcga_tpm_log2,
     tme1.list1,
     tme1.list2,
     tme2.list,
     ImmuCellAI.list,
     ssGSEAScore_by_muti_group_genes,
     file = 'ssgsea.Rdata')
load('ssgsea.Rdata')

tcga_tme1.list1_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga_tpm_log2,
                                                         genelist = tme1.list1)
tcga_tme1.list2_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga_tpm_log2,
                                                         genelist = tme1.list2)
tcga_tme2.list_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga_tpm_log2,
                                                        genelist = tme2.list)
tcga_ImmuCellAI.list_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga_tpm_log2,
                                                              genelist = ImmuCellAI.list)



save(tcga_tme1.list1_score,
     tcga_tme1.list2_score,
     tcga_tme2.list_score,
     tcga_ImmuCellAI.list_score,
     file = 'ssgsea_res.Rdata')
load('ssgsea_res.Rdata')

tcga_tme1.list1_score["Activated B cell", rownames(tcga_model_data)]

tcga_tme1.list1_score_plot <- wb_boxplot(dat = tcga_tme1.list1_score[, rownames(tcga_model_data)],
                                         groups = tcga_model_data$RiskType,
                                         title = 'RiskType',
                                         xangle = 60,
                                         ylab = 'Score',
                                         col = mycolor)
tcga_tme1.list1_score_plot          # 
# 
tcga_tme1.list1_score <- as.data.frame(t(tcga_tme1.list1_score[, rownames(tcga_model_data)]))
tcga_hs_RS_cor <- psych::corr.test(as.matrix(tcga_tme1.list1_score))
tcga_hs_RS_cor_R <- tcga_hs_RS_cor$r
tcga_hs_RS_cor_P <- tcga_hs_RS_cor$p
library(ggcorrplot)
tcga_hs_RS_cor_plot <- ggcorrplot(tcga_hs_RS_cor_R, 
                                  hc.order = TRUE,
                                  colors = c(mycolor[1], "white", mycolor[2]),
                                  ggtheme = ggplot2::theme_bw,
                                  # type = c("lower"),
                                  p.mat = tcga_hs_RS_cor_P)
tcga_hs_RS_cor_plot
ggsave(plot = tcga_hs_RS_cor_plot,
       filename = '05_RS_immune/tcga_hs_RS_cor_plot.pdf',
       width = 8, height = 8)

tcga_tme1.list2_score_plot <- wb_boxplot(dat = tcga_tme1.list2_score[, rownames(tcga_model_data)],
                                         groups = tcga_model_data$RiskType,
                                         xangle = 0,
                                         title = 'RiskType',
                                         ylab = 'Score',
                                         col = mycolor)
tcga_tme1.list2_score_plot           # 
tcga_tme1.list2_score_list <- list()
for (fea in rownames(tcga_tme1.list2_score)) {
  print(fea)
  tmp <- mg_violin(data.frame(tcga_model_data$RiskType,
                              as.numeric(tcga_tme1.list2_score[fea, 
                                                               rownames(tcga_model_data)])), 
                   melt=TRUE, ylab=fea,
                   leg.title='', 
                   jitter = F,
                   xlab = 'RiskType',
                   legend.pos='tl',
                   show_compare = T)
  tcga_tme1.list2_score_list[[fea]] <- tmp
}

tcga_tme1.list2_score_plot <- cowplot::plot_grid(plotlist = tcga_tme1.list2_score_list,
                                                 ncol = 2)
tcga_tme1.list2_score_plot     #



tcga_estimate <- immu_estimate(exp = tcga_tpm_log2)
tcga_estimate_clister_list <- list()
for (fea in colnames(tcga_estimate)) {
  print(fea)
  tmp <- mg_violin(data.frame(tcga_model_data$RiskType,
                              as.numeric(tcga_estimate[rownames(tcga_model_data), fea])), 
                   melt=TRUE, ylab=fea,
                   leg.title='', 
                   jitter = F,
                   xlab = 'RiskType',
                   legend.pos='tl',
                   show_compare = T)
  tcga_estimate_clister_list[[fea]] <- tmp
}

tcga_estimate_cluster_plot <- cowplot::plot_grid(plotlist = tcga_estimate_clister_list,
                                                 ncol = 3)
tcga_estimate_cluster_plot          #

ggsave(plot = tcga_estimate_cluster_plot,
       filename = '05_RS_immune/tcga_estimate_cluster_plot.pdf',
       width = 9, height = 5)


tcga_mcpcounter <- immu_MCPcounter(exp = tcga_tpm_log2)
dim(tcga_mcpcounter)

tcga_mcpcounter_score_plot <- wb_boxplot(dat = t(tcga_mcpcounter[rownames(tcga_model_data), ]),
                                         groups = tcga_model_data$RiskType,
                                         title = 'RiskType',
                                         xangle = 60,
                                         ylab = 'Score',
                                         col = mycolor)


tcga_mcpcounter_score_plot    #


tcga_tme2.list_score_plot <- wb_boxplot(dat = tcga_tme2.list_score[, rownames(tcga_model_data)],
                                        groups = tcga_model_data$RiskType,
                                        title = 'RiskType',
                                        xangle = 60,
                                        ylab = 'Score',
                                        col = mycolor)
tcga_tme2.list_score_plot

tcga_ImmuCellAI.list_score_plot <- wb_boxplot(dat = tcga_ImmuCellAI.list_score[, rownames(tcga_model_data)],
                                              groups = tcga_model_data$RiskType,
                                              title = 'RiskType',
                                              xangle = 60,
                                              ylab = 'Score',
                                              col = mycolor)
tcga_ImmuCellAI.list_score_plot

tcga_immune_cluster_plot1 <- cowplot::plot_grid(tcga_hs_RS_cor_plot,
                                                tcga_tme1.list2_score_plot,
                                                ncol = 2,
                                                labels = LETTERS[2:3],
                                                rel_widths = c(3,2))
tcga_immune_cluster_plot1

tcga_immune_cluster_plot2 <- cowplot::plot_grid(tcga_estimate_cluster_plot,
                                                tcga_mcpcounter_score_plot,
                                                ncol = 2,
                                                labels = LETTERS[4:5],
                                                rel_widths = c(3,2))
tcga_immune_cluster_plot2

tcga_immune_cluster_plot <- cowplot::plot_grid(tcga_tme1.list1_score_plot,
                                               tcga_immune_cluster_plot1,
                                               tcga_immune_cluster_plot2,
                                               ncol = 1,
                                               labels = 'A',
                                               label_fontfamily = 'Times',
                                               label_size = 14)
tcga_immune_cluster_plot
ggsave(plot = tcga_immune_cluster_plot,
       filename = '05_RS_immune/fig5.pdf',
       width = 18, height = 18)


#  ###################
dir.create('06_TIDE')
tcga_tide_dat <- t(scale(t(tcga_tpm_log2[, rownames(tcga_model_data)]),scale = F))

write.table(tcga_tide_dat,
            file = '06_TIDE/tcga_tide_dat.txt',
            quote = F, sep = '\t')

tcga_tide_dat<-read.csv('06_TIDE/tcga_tide_res.csv',row.names = 1,stringsAsFactors = F)
colnames(tcga_tide_dat)


colnames(tcga_tide_dat)
tcga_tide_RSplot1 <- wb_boxplot(dat = t(tcga_tide_dat[rownames(tcga_model_data), c("TIDE",
                                                                                 "Dysfunction",
                                                                                 "Exclusion")]),
                               groups = tcga_model_data$RiskType,
                               title = 'RiskType',
                               xangle = 0,
                               ylab = 'Score',
                               col = mycolor)
tcga_tide_RSplot1

tcga_tide_list <- list()
for (fea in c("TIDE",
              "Dysfunction",
              "Exclusion")) {
  print(fea)
  tmp <- mg_violin(data.frame(tcga_model_data$RiskType,
                              as.numeric(tcga_tide_dat[rownames(tcga_model_data), fea])), 
                   melt=TRUE, ylab=fea,
                   leg.title='', 
                   jitter = F,
                   xlab = 'RiskType',
                   legend.pos='tl',
                   show_compare = T)
  tcga_tide_list[[fea]] <- tmp
}

tcga_tide_RSplot <- cowplot::plot_grid(plotlist = tcga_tide_list,
                                       ncol = 3,
                                       labels = LETTERS[1:3])
tcga_tide_RSplot
# ggsave(plot = tcga_tide_RSplot,
#        filename = '06_TIDE/tcga_tide_RSplot.pdf',
#        width = 9, height = 5)



tcga_tide_dat <- cbind(tcga_tide_dat,
                       tcga_model_data[rownames(tcga_tide_dat), c("RS", "RiskType")])


tcga_tide_rs_cor <- cor_point(x = tcga_tide_dat$TIDE,
                              y = tcga_tide_dat$RS,
                              xlab = 'TIDE',
                              ylab = 'RiskScore',
                              top_col = mycolor[1],
                              right_col = mycolor[2])
tcga_tide_rs_cor

tcga_Dysfunction_rs_cor <- cor_point(x = tcga_tide_dat$Dysfunction,
                                     y = tcga_tide_dat$RS,
                                     xlab = 'Dysfunction',
                                     ylab = 'RiskScore',
                                     top_col = mycolor[1],
                                     right_col = mycolor[2])
tcga_Dysfunction_rs_cor

tcga_Exclusion_rs_cor <- cor_point(x = tcga_tide_dat$Exclusion,
                                   y = tcga_tide_dat$RS,
                                   xlab = 'Exclusion',
                                   ylab = 'RiskScore',
                                   top_col = mycolor[1],
                                   right_col = mycolor[2])
tcga_Exclusion_rs_cor

tcga_tide_plot <- cowplot::plot_grid(tcga_tide_rs_cor,
                                     tcga_Dysfunction_rs_cor,
                                     tcga_Exclusion_rs_cor,
                                     labels = LETTERS[4:6],
                                     ncol = 3)
tcga_tide_plot

fig6 <- cowplot::plot_grid(tcga_tide_RSplot,
                           tcga_tide_plot,
                           ncol = 1,
                           label_fontfamily = 'Times',
                           label_size = 14)
fig6

ggsave(plot = fig6,
       filename = '06_TIDE/fig6.pdf',
       width = 12, height = 8)


#  ####
dir.create('07_drugs')
##  
# 
gdsc_exp <- read.delim('00_origin_datas/GDSC/Cell_line_RMA_proc_basalExp.xlsx',header = T, check.names = F)
dim(gdsc_exp)   # 17737 1020

gdsc_exp <- gdsc_exp[which(gdsc_exp$GENE_SYMBOLS != ''),]
dim(gdsc_exp)  # 17419 1020
rownames(gdsc_exp) <- gdsc_exp$GENE_SYMBOLS
gdsc_exp <- gdsc_exp[,c(-1,-2)]
range(gdsc_exp)

## 
gdsc_infor <- read.csv('00_origin_datas/GDSC/LUAD_IC_Tue Jun  6 06_55_41 2023.csv')
dim(gdsc_infor)    #  15747 13
gdsc_infor$Cosmic.ID <- paste('DATA.',gdsc_infor$Cosmic.ID,sep = '')


gongyou_gdsc <- intersect(colnames(gdsc_exp),gdsc_infor$Cosmic.ID)
gongyou_gdsc

gdsc_exp <- gdsc_exp[,gongyou_gdsc]

gdsc.genes <- intersect(genes, rownames(gdsc_exp))
gdsc.genes
genes
gdsc_lan <- lan
gdsc_lan

risk.cellline <- as.matrix(t(gdsc_lan%*%as.matrix(gdsc_exp[gdsc.genes,])))
dim(risk.cellline)    # 61 1

risk.cellline<-data.frame(RiskScore = risk.cellline)
risk.cellline$samples <- rownames(risk.cellline)

colnames(gdsc_infor)
gdsc_drug_ic_auc <- gdsc_infor[,c('Cosmic.ID','Drug.Name','IC50','AUC')]
colnames(gdsc_drug_ic_auc)[1] <- 'samples'

gdsc.risk_drug<-merge(gdsc_drug_ic_auc,risk.cellline,by='samples')
colnames(gdsc.risk_drug)[5] <- 'RiskScore'

dim(gdsc.risk_drug)   # 15567  5
head(gdsc.risk_drug)
gdsc.risk_drug=crbind2DataFrame(gdsc.risk_drug)

## 
library(dplyr)
gdsc.risk_drug.cor=gdsc.risk_drug %>% 
  dplyr::select(Drug.Name, AUC, RiskScore) %>%
  dplyr::group_by(Drug.Name) %>% 
  dplyr::summarize(cor=psych::corr.test(AUC,RiskScore,method = "spearman",adjust = "none")$r,pval=psych::corr.test(AUC,RiskScore,method = "spearman",adjust = "none")$p)
dim(gdsc.risk_drug.cor)    # 288 3

gdsc.risk_drug.cor$Class=case_when(gdsc.risk_drug.cor$cor>0.3 & gdsc.risk_drug.cor$pval<0.05 ~ "positive",
                                   gdsc.risk_drug.cor$cor < (-0.3) & gdsc.risk_drug.cor$pval < 0.05 ~ "negative",
                                   TRUE ~ "ns" )
gdsc.risk_drug.cor=data.frame(gdsc.risk_drug.cor)
table(gdsc.risk_drug.cor$Class)
head(gdsc.risk_drug.cor)

########## 
gdsc_drug_selected<-gdsc.risk_drug.cor[which(gdsc.risk_drug.cor$Class != 'ns'),]
gdsc_drug_selected<-gdsc_drug_selected[order(gdsc_drug_selected$cor),]
gdsc_drug_selected$`-log10(pvalue)`<-(-log10(gdsc_drug_selected$pval))
head(gdsc_drug_selected)

library(tidyverse)
fig7A=gdsc_drug_selected %>% 
  ggplot(aes(reorder(Drug.Name, cor), cor)) + 
  geom_col(aes(fill = `-log10(pvalue)`)) + 
  scale_fill_gradient2(low = mycolor[1],
                       high = mycolor[2],
                       midpoint = 1.3) +
  #coord_flip() + 
  labs(x = "")+
  labs(y = "Rs of drug sensitivity and RiskScore")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
fig7A
savePDF('07_drugs/fig7A.pdf',fig7A,height = 5,width = 6)



####  
head(gdsc.risk_drug)
gdsc_drug_selected$Drug.Name
# "GSK2256098C" "Carmustine"  "Dacarbazine" "PCI-34051"   "Elephantin" 
gdsc.risk_drug_dat <- data.frame()
for (drug in gdsc_drug_selected$Drug.Name) {
  print(drug)
  tmp <- gdsc.risk_drug[which(gdsc.risk_drug$Drug.Name == drug),]
  tmp$RiskGroup <- ifelse(tmp$RiskScore>median(tmp$RiskScore),'High','Low')
  gdsc.risk_drug_dat <- rbind(gdsc.risk_drug_dat, tmp)
}
head(gdsc.risk_drug_dat)
table(gdsc.risk_drug_dat$DRUG_NAME)

library(ggpubr)
gdsc.risk_drug_plot <- ggplot(gdsc.risk_drug_dat, 
                              aes(x=Drug.Name, y=AUC, fill=RiskGroup)) +
  geom_boxplot(notch = T) +  
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic() +scale_fill_manual(values = mycolor)+
  theme(axis.text.x = element_text(angle=90, hjust = 1,vjust = 0.5)) +
  xlab('Drugs') + ylab('AUC')
gdsc.risk_drug_plot
ggsave(plot = gdsc.risk_drug_plot,
       filename = '07_drugs/fig7B.pdf',
       width = 6, height = 5)

fig7 <- cowplot::plot_grid(fig7A,
                        gdsc.risk_drug_plot,
                        ncol = 2,
                        labels=c('A','B'),
                        label_fontfamily = 'Times',
                        label_size = 14)
fig7
ggsave(plot = fig7,
       filename = '07_drugs/fig7.pdf',
       width = 12, height = 5)

#  ##############
dir.create('08_nomogram')
library(rpart)
library(rpart.plot)

all.t.cli.forModel=tcga_RS_cli
colnames(all.t.cli.forModel)

table(all.t.cli.forModel$Stage)
all.t.cli.forModel$Stage=ifelse(all.t.cli.forModel$Stage %in% c('I','II'),'I+II','III+IV')
all.t.cli.forModel$Stage=as.factor(all.t.cli.forModel$Stage)

table(all.t.cli.forModel$T.Stage)
all.t.cli.forModel$T.Stage=ifelse(all.t.cli.forModel$T.Stage %in% c('T1','T2'),'T1+T2','T3+T4')
all.t.cli.forModel$T.Stage=as.factor(all.t.cli.forModel$T.Stage)


table(all.t.cli.forModel$N.Stage)
all.t.cli.forModel$N.Stage=ifelse(all.t.cli.forModel$N.Stage %in% c('N0','N1'),'N0+N1','N2+N3')
all.t.cli.forModel$N.Stage=as.factor(all.t.cli.forModel$N.Stage)

table(all.t.cli.forModel$M.Stage)
all.t.cli.forModel$M.Stage=as.factor(all.t.cli.forModel$M.Stage)


table(all.t.cli.forModel$Age)
all.t.cli.forModel$Age <- ifelse(all.t.cli.forModel$Age > 65, '>65', '<=65')
all.t.cli.forModel$Age=factor(all.t.cli.forModel$Age,levels = c('<=65','>65'))

table(all.t.cli.forModel$Gender)


str(all.t.cli.forModel)
colnames(all.t.cli.forModel)
library(survival)
library(ggsci)
pfit <- rpart(Surv(OS.time, OS) ~ Age + Gender + T.Stage + N.Stage + M.Stage + Stage + RiskType, 
              data = all.t.cli.forModel)
print(pfit)

printcp(pfit)

pfit2 <- prune(pfit, cp=0.015)
print(pfit2)
rpart.plot(pfit2)


pdf('08_nomogram/fig8A.pdf',height = 5,width = 4)
rpart.plot(pfit2)
dev.off()

all.t.cli.forModel$class=pfit2$where
table(all.t.cli.forModel$class)

all.t.cli.forModel$class[all.t.cli.forModel$class==3]='C1'
all.t.cli.forModel$class[all.t.cli.forModel$class==4]='C2'
all.t.cli.forModel$class[all.t.cli.forModel$class==5]='C3'


table(all.t.cli.forModel$RiskType,all.t.cli.forModel$class)

fit <- survfit( Surv(OS.time/365, OS) ~ class,data = all.t.cli.forModel)
library(survminer)
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
fig8b=ggsurvplot(fit,data=all.t.cli.forModel,
                  conf.int = T,
                  pval = TRUE,
                  fun = "pct",
                  risk.table = TRUE,
                  size = 1,
                  title='',
                  ggtheme=theme_bw(),
                  # linetype = "strata",
                  palette = pal_d3()(9)[1:5],
                  legend = "bottom",
                  legend.title = "Cluster",
                  legend.labs = c("C1","C2",'C3')
)
fig8b
fig8B=ggpubr::ggarrange(fig8b$plot,fig8b$table, ncol = 1, nrow = 2
                         ,heights = c(0.7,0.3)
                         ,align = "v",common.legend = T)
fig8B
table(all.t.cli.forModel$RiskType
      ,all.t.cli.forModel$class)

fig8c=plotMutiBar(table(all.t.cli.forModel$RiskType
                         ,all.t.cli.forModel$class))
fig8c
fig8d=plotMutiBar(table(all.t.cli.forModel$OS
                         ,all.t.cli.forModel$class))
fig8d

fig8bcd=mg_merge_plot(fig8B,fig8c,fig8d,nrow = 1,ncol = 3,labels = c('B','C','D'))
fig8bcd

savePDF('08_nomogram/fig8BCD.pdf',fig8bcd,height = 6,width = 12)

# ########
mg_compare_uni_muti_cox_use=function(dat,event,os){
  # dat=crbind2DataFrame(dat)
  sig.clini=colnames(dat)
  dat$time=os
  dat$status=event
  dat=dat[which(!is.na(os)&!is.na(event)),]
  all.cox=rbind()
  rnames=c()
  for(s in sig.clini){
    fmla <- as.formula(paste0("Surv(time, status) ~",s))
    cox <- coxph(fmla, data = dat)
    #summary(cox)[[7]]
    #print(summary(cox))
    re=cbind(summary(cox)[[7]][,5],summary(cox)[[7]][,2],summary(cox)[[8]][,3],summary(cox)[[8]][,4])
    if(nrow(re)==1){
      rnames=c(rnames,s)
    }else{
      rnames=c(rnames,row.names(summary(cox)[[7]]))
    }
    all.cox=rbind(all.cox,re)
  }
  row.names(all.cox)=rnames
  colnames(all.cox)=c('p.value','HR','Low 95%CI','High 95%CI')
  
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(sig.clini,collapse = '+')))
  cox <- coxph(fmla, data = dat)
  muti.re=cbind(summary(cox)[[7]][,5],summary(cox)[[7]][,2],summary(cox)[[8]][,3],summary(cox)[[8]][,4])
  row.names(muti.re)=row.names(summary(cox)[[7]])
  colnames(muti.re)=c('p.value','HR','Low 95%CI','High 95%CI')
  return(list(muti=crbind2DataFrame(muti.re),uni=crbind2DataFrame(all.cox)))
}
getForestplotData=function(res_sig){
  res_sig<-signif(res_sig,digits=2)
  res_sig$CI_for_HR=paste0(" (",res_sig$`Low 95%CI`, "-", res_sig$`High 95%CI`, ")")
  colnames(res_sig)=c("p.value","HR","CI_lower","CI_upper", "(95%_CI_for_HR)")
  
  forest_table<-data.frame(Features=rownames(res_sig),HR=res_sig$HR,`(95%CI)`=res_sig$`(95%_CI_for_HR)`, pvalue=res_sig$p.value,check.names = F,stringsAsFactors = F)
  forest_table$sig<-mg_format_p_values(forest_table$pvalue)
  
  forest_table2<-data.frame(Features="Features",HR="HR",`(95%CI)`="(95%CI)", pvalue="p-value",sig='Significant',check.names = F,stringsAsFactors=F)
  tabletext<-rbind(forest_table2,forest_table)
  
  forest_stastic<-data.frame(mean=as.numeric(as.character(res_sig$HR)),lower=as.numeric(as.character(res_sig$CI_lower)),upper=as.numeric(as.character(res_sig$CI_upper)))
  forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
  cochrane_from_rmeta<-rbind(forest_stastic1,forest_stastic)
  return(list(tabletext,cochrane_from_rmeta))
}

library(forestplot)

all.t.cli.forCox=all.t.cli.forModel
colnames(all.t.cli.forCox)
colnames(all.t.cli.forCox)[c(17,11:16)]

cm.cox=mg_compare_uni_muti_cox_use(all.t.cli.forCox[,c(17,11:16)],
                                   os=all.t.cli.forCox$OS.time,
                                   event = all.t.cli.forCox$OS)

cm.cox.uni=signif(cm.cox$uni,digits=3)
cm.cox.uni$`Hazard Ratio(95%CI)`=paste0(cm.cox.uni$HR,"(",cm.cox.uni$`Low 95%CI`, "-", cm.cox.uni$`High 95%CI`, ")")

table(cm.cox.uni$p.value<0.05)
cm.cox.uni[which(cm.cox.uni$p.value<0.05),]

####### 
colnames(all.t.cli.forCox)[c(17,13:16)]
cm.cox2=mg_compare_uni_muti_cox_use(all.t.cli.forCox[,c(17,13:16)],
                                    os=all.t.cli.forCox$OS.time,
                                    event = all.t.cli.forCox$OS)
cm.cox.muti=signif(cm.cox2$muti,digits=3)
cm.cox.muti
cm.cox.muti$`Hazard Ratio(95%CI)`=paste0(cm.cox.muti$HR,"(",cm.cox.muti$`Low 95%CI`, "-", cm.cox.muti$`High 95%CI`, ")")

cm.cox.uni=getForestplotData(cm.cox$uni)
cm.cox.muti=getForestplotData(cm.cox2$muti)

pdf('08_nomogram/fig8E.pdf',height = 6,width = 8,onefile = F)
tabletext=cm.cox.uni[[1]]
cochrane_from_rmeta=cm.cox.uni[[2]]
forestplot(tabletext, 
           mean=cochrane_from_rmeta$mean,
           lower=cochrane_from_rmeta$lower,
           upper=cochrane_from_rmeta$upper,
           graph.pos=2,
           hrzl_lines=list('2'=gpar(lty=1,col="black"),
                           '9'=gpar(lty=1,col="black")),
           zero = 1, 
           xlog=T,
           fn.ci_norm = fpDrawDiamondCI,
           boxsize = 0.2, 
           col=fpColors(line = mycolor[1], 
                        box=mycolor[2], 
                        zero = mycolor[3], 
                        summary=mycolor[4]), 
           lty.ci = 7,  
           lwd.ci = 1,
           ci.vertices.height = 0.05, 
           lineheight = "auto", 
           xlab="Hazard ratio" )
dev.off()

pdf('08_nomogram/fig8F.pdf',height = 6,width = 8,onefile = F)
tabletext=cm.cox.muti[[1]]
cochrane_from_rmeta=cm.cox.muti[[2]]
forestplot(tabletext, 
           mean=cochrane_from_rmeta$mean,
           lower=cochrane_from_rmeta$lower,
           upper=cochrane_from_rmeta$upper,
           graph.pos=2,
           hrzl_lines=list('2'=gpar(lty=1,col="black"),
                           '7'=gpar(lty=1,col="black")),
           zero = 1, 
           xlog=T,
           fn.ci_norm = fpDrawDiamondCI,
           boxsize = 0.2, 
           col=fpColors(line = mycolor[1], 
                        box=mycolor[2], 
                        zero = mycolor[3], 
                        summary=mycolor[4]), 
           lty.ci = 3,  
           lwd.ci = 1,
           ci.vertices.height = 0.05, 
           lineheight = "auto", 
           xlab="Hazard ratio" )

dev.off()

mg_plotDCA=function(status,fmlas,modelNames,data){
  set.seed(123)
  all.mod=list()
  for(i in 1:length(fmlas)){
    fmla <- as.formula(paste0("status~",fmlas[i]))
    model<-rmda::decision_curve(fmla,
                                data=data,
                                bootstraps=500)
    all.mod=c(all.mod,list(model))
  }
  rmda::plot_decision_curve(all.mod,
                            curve.names=modelNames,
                            # col=mg_colors[c(1,10:12,4,5,7,8)],
                            col=mg_colors[c(1,4,5,8)],
                            xlim=c(0,1),legend.position="topright",
                            lwd=1,
                            confidence.intervals=FALSE)
}

colnames(all.t.cli.forCox)[c(17,13)]


mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',quick=T){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=c(1,3,5)
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*1,12*3,12*5)
  }else{
    cut.time=c(365*1,365*3,5*365)
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['95%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#
  #,observation=pbc[2,] #
  #
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #
  #              ,showP = T #
  #              ,droplines = F#
  #,colors = mg_colors[1:3] #
  #,rank="decreasing") #
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c('1-Year Survival','3-Year survival'
                           ,'5-Year survival')[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c("1-year","3-year",'5-year')[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}


pdf('08_nomogram/nomogram.pdf', width = 12 , height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore = all.t.cli.forCox$RiskScore,
                                T.Stage = all.t.cli.forCox$T.Stage),
                     os = all.t.cli.forCox$OS.time,
                     status = all.t.cli.forCox$OS)
dev.off()



# #####
save.image('LUAD_001.Rdata')
load('LUAD_001.Rdata')












