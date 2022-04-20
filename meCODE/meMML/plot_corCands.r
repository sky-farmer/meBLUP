################################################################################
#################################plot_corCands.r################################
################################################################################
# R-Script to create Polts for output
#
#-------------------------------------------
#USAGE:
#       plot_corCan(EBVcan)
#
#
#--------------------------------------------
#INPUT:
#   EBVcan.....data.table with colunms ID and EBVs
#
#-----------------------------------------------------------------
#OUTPUT:
#   CorPlot......
#   
#
#
# Himmelbauer, 3.1.2022
##################################################################################

plot_corCan=function(EBVcan, x_hist, RESID){

    #--------------------------------------------
    # Plot Correaltions between Solutions and TBVs
    #--------------------------------------------
    tbv<-          fread(paste0(PATH,'.bv'))
    print(paste0(  'Iterationen: ',(ncol(EBVcan)-1)))
    setkey(EBVcan, 'ID')
    setkey(tbv,    'id')
    tbv<-          tbv[J(EBVcan$ID)]
    correl<-       NULL
    
    for (i in (2:(ncol(EBVcan)))) {
        c<-        cor(tbv$true, EBVcan[,..i])
        correl<-   c(correl, c)
    }
    correl<-       data.table(Iteration=1:length(correl), Correlation=correl, Rep=x_hist$Rep)
    
    #Plot correlations EBVs vs. TBVs for test animals
    p <-           ggplot(correl, aes(Iteration, Correlation))+
                        geom_line(size=0.8)+
                        ylim(0,1)+ 
                        theme(axis.text=element_text(size=14),
                              axis.title=element_text(size=16,face="bold"))
    png(paste0(PATH,'CorPlot_Cands',suffix,'.png'),width=12*300,height=7*300,res=300);plot(p);dev.off()
    if (PRT==1) {
        fwrite(correl,   paste0(PATH,'.corr_test',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    }
    
    print(paste0('highest Correlation TBV-meBLUP at iteration:', correl[correl$Correlation==max(correl$Correlation),c('Iteration')]))
    print(paste0('highest Correlation TBV-meBLUP:',              correl[correl$Correlation==max(correl$Correlation),c('Correlation')]))
    
    ITB1<<-        correl$Iteration[correl$Correlation==max(correl$Correlation)][1]
    BEST<-         ITB1+1
    best<<-        EBVcan[,..BEST]
    correl2<-      NULL
    for (i in (2:(ncol(EBVcan)))) {
        c<-        cor(best, EBVcan[,..i])
        correl2<-  c(correl2, c)
    }
    correl2<-      data.table(Iteration=1:length(correl2), Correlation=correl2, Rep=x_hist$Rep)
    ITB4<<-        correl2$Iteration[correl2$Correlation>=0.99][1]
    
    correlsst<-    NULL
    for (i in (1:(ncol(EBVcan)-1))) {
        c<-        cor(tbv$sstep, EBVcan[,..i])
        correlsst<-c(correlsst, c)
    }
    correlsst<-    data.table(Iteration=1:length(correlsst), Correlation=correlsst, Rep=x_hist$Rep)
    
    #Plot correlations EBVs vs. EBVs from sstep for test animals
    p <-           ggplot(correlsst, aes(Iteration, Correlation))+
                        geom_line(size=0.8)+
                        ylim(0,1)+ 
                        theme(axis.text=element_text(size=14),
                              axis.title=element_text(size=16,face="bold"))
    png(paste0(PATH,'CorPlot_sstep_Cands',suffix,'.png'),width=12*300,height=7*300,res=300);plot(p);dev.off()
    
    print(paste0('highest Correlation sstep-meBLUP at iteration:', correlsst[correlsst$Correlation==max(correlsst$Correlation),c('Iteration')]))
    print(paste0('highest Correlation sstep-meBLUP:', correlsst[correlsst$Correlation==max(correlsst$Correlation),c('Correlation')]))
    print(paste0('Correlation TBV-sstep:', cor(tbv$true, tbv$sstep)))

    if (RESID) {
        par<-           as.matrix(x_hist)
        RES<-           NULL
        for (i in (1:(nrow(x_hist)))) {
            xi<-        par[i,-1]
            u<-         xi[(nf+1):(nf+na)]
            r<-         xi[(nf+na+nSNP+1)]
            i<-         i+1
            ucan<-      as.vector(EBVcan[,..i]$EBV)
            Rdcan<-     data.table(ID=A$can, resid=rep(1, length(A$can)))
            Rdcan$resid[Rdcan$ID%in%A$can]<-r
            Rdcan<-     Rdcan$resid
            uPA<-       A$PAcan%*%u
            res<-       sd(ucan-uPA)
            RES<-       c(RES, res)
        }
  
        d<-             data.table(Iteration=correl$Iteration , Correlation=correl$Correlation, logRES=log10(RES))
        coeff<-         max(d$logRES,na.rm=TRUE)*1.1
    
        print(paste0('Correlation TBV-meBLUP at iteration with best log(sd(u-PA)):', d[d$logRES==min(d$logRES),c('Correlation')]))
        print(paste0('Iteration with best log(sd(u-PA)):', d[d$logRES==min(d$logRES),c('Iteration')]))
        ITB2<<-        d$Iteration[d$logRES==min(d$logRES,na.rm=TRUE)][1]
    
        #Plot Resid and correaltion between EBVs and TBVs in one plot
        p<-            ggplot(d, aes(x=Iteration)) +
                            geom_line( aes(y=Correlation, linetype='Correlation'),    size=0.8) + 
                            geom_line( aes(y=logRES/coeff, linetype='log(sd(u-PA))'), size=0.8, ) +
                            scale_linetype_manual(values=c('solid','dashed'),name='') +
                            geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                            geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                            geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                            scale_color_manual(name="", values=c("red","green","blue"))+
                            scale_y_continuous(
                                name = "Correlation (Test)",
                                sec.axis = sec_axis(~.*coeff, name="log(sd(u-PA))")) + 
                            theme(axis.text=element_text(size=14),
                                  axis.title=element_text(size=16,face="bold"),
                                  legend.title = element_text(size=16), #change legend title font size
                                  legend.text = element_text(size=16),
                                  legend.position="bottom",
                                  legend.key.width = unit(1.9, "cm"),
                                  legend.key.height = unit(0.7, "cm"),
                                  legend.box="vertical", legend.margin=margin())+
                                  guides(color     = guide_legend(order = 0),
                                         linetype  = guide_legend(order = 1))
        png(paste0(PATH,'RESIDCorPlot_Cands',suffix,'.png'),width=12*300,height=7*300,res=300);plot(p);dev.off()
        if (PRT==1) {
            fwrite(d,   paste0(PATH,'.RESIDcorr_test',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
        }
    }
}




