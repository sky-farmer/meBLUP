################################################################################
#################################plot_out.r#####################################
################################################################################
# R-Script to create Polts for output
#
#-------------------------------------------
#USAGE:
#       plot_out(x_hist, PATH, suffix)
#
#
#--------------------------------------------
#INPUT:
#   x.............historic solution vectors
#   PATH..........full path to files
#   suffix........added to the end of the filenames
#
#-----------------------------------------------------------------
#OUTPUT:
#   CorPlot......
#   h2Plot.......
#   sigsnpPlot...
#   rPlot........ 
#
#
# Himmelbauer, 28.12.2021
##################################################################################
library(data.table)
library(ggplot2)
library(Matrix)
library(runner)

plot_out=function(x_hist, PATH, SUFFIX){

    #--------------------------------------------
    # Plot Correaltions between Solutions and TBVs
    #--------------------------------------------
    tbv<-          fread(paste0(PATH,'.bv'))
    ped<-          fread(paste0(PATH,'.ped'), fill=TRUE, header=FALSE)
    ped<-          ped[ped$V1%in%A$cal]
    tbv<-          tbv[tbv$id%in%A$cal]
    print(paste0(  'Iterationen: ',nrow(x_hist)))
    ebv<-          as.data.table(cbind(ped[,1], t(x_hist[,(1+nf+1):(1+nf+na)])))
    names(ebv)[1]<-'ID'
    setkey(ebv,    'ID')
    setkey(tbv,    'id')
    ebv<-          ebv[J(tbv$id)]
    correl<-       NULL
    
    for (i in (2:(ncol(ebv)))) {
        c<-cor(tbv$true, ebv[,..i])
        correl<-   c(correl, c)
    }
    correl<-       data.table(Iteration=1:length(correl), Correlation=correl, Rep=x_hist$Rep)
    
    print(paste0('highest Correlation_training TBV-meBLUP at iteration:', correl[correl$Correlation==max(correl$Correlation),c('Iteration')]))
    print(paste0('highest Correlation_training TBV-meBLUP:', correl[correl$Correlation==max(correl$Correlation),c('Correlation')]))
    print(paste0('Correlation TBV-sstep:', cor(tbv$true, tbv$sstep)))
    
    #Plot correlations EBVs vs. TBVs for training animals
    p <-            ggplot(correl, aes(Iteration, Correlation))+
                        geom_line(size=0.8)+
                        ylim(0,1)+
                        geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                        geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                        scale_color_manual(name="", values=c("red","green"))+
                        ylab('Correlation (Training)')+
                        theme(axis.text=element_text(size=14), 
                              axis.title=element_text(size=16,face="bold"),
                              legend.title = element_text(size=16),
                                      legend.text = element_text(size=16),
                                      legend.position="bottom",
                                      legend.key.width = unit(1.9, "cm"),
                                      legend.key.height = unit(0.7, "cm"),
                                      legend.box="vertical", legend.margin=margin())+
    png(paste0(PATH,'CorPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(p);dev.off()
    if (PRT==1) {
        fwrite(correl,   paste0(PATH,'.corr_train',SUFFIX), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    }
    
    
    if (PRT==1) {
    if ('estim_r' %in% names(x_hist)) {
        #--------------------------------------------
        # Plot estim_r for all iterations
        #--------------------------------------------
        erhist<-           x_hist[,c('r','estim_r')]
        erhist$Iteration<- 1:nrow(erhist)
        d <-               melt(erhist, id.vars="Iteration")
        
        erp <-             ggplot(d,aes(Iteration, value, linetype=variable))+
                                geom_line(size=0.8)+
                                scale_linetype_manual(values=c('solid','dotted'), labels=c(expression(paste(rho)),expression(paste(rho, '_estim'))),name='')+
                                scale_size_manual(values=c(0.8,0.8), labels=c(expression(paste(rho)),expression(paste(rho, '_estim'))),name='')+
                                #ylim(0,0.6)+
                                ylim(0,1.1*max(d$value[d$Iteration>20]))+
                                ylab(expression(paste(rho)))+
                                #geom_text(aes(max(d$Iteration),0.38,label = 'optim=0.38', vjust = -1, hjust=+1.0))+
                                #geom_hline(yintercept=0.38, linetype='dashed', color='grey25', size=1.0)+
                                geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                scale_color_manual(name="", values=c("red","green","blue"))+
                                theme(axis.text=element_text(size=14),
                                      axis.title=element_text(size=16,face="bold"),
                                      legend.title = element_text(size=16),
                                      legend.text = element_text(size=16),
                                      legend.position="bottom",
                                      legend.key.width = unit(1.9, "cm"),
                                      legend.key.height = unit(0.7, "cm"),
                                      legend.box="vertical", legend.margin=margin())+
                                guides(color     = guide_legend(order = 0),
                                       linetype  = guide_legend(order = 1),
                                       size      = guide_legend(order = 1))
        png(paste0(PATH,'er_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(erp);dev.off()
        
        #--------------------------------------------
        # Plot estim_h2_1/2/3 for all iterations
        #--------------------------------------------
        eh2hist<-           x_hist[,c('h2', 'h21', 'h22', 'h23', 'h24', 'h21c', 'h22c', 'h23c', 'h24c')]
        eh2hist$h2<-        eh2hist$h2^2
        eh2hist$Iteration<- 1:nrow(eh2hist)
        d <-                melt(eh2hist, id.vars="Iteration")
        
        eh2p <-             ggplot(d, aes(Iteration, log10(value), linetype=variable)) +
                                geom_line(aes(size=variable))+
                                scale_linetype_manual(values=c('solid','dotted','dashed','dotdash', 'twodash','dotted','dashed','dotdash', 'twodash'),
                                    labels=c('h2', 'h21(PHE)', 'h22(PED)', 'h23(HYB)','h24(PED+HYB)', 'c*h21(PHE)', 'c*h22(PED)', 'c*h23(HYB)', 'c*h24(PED+HYB)'),
                                    name='')+
                                scale_size_manual(values=c(0.8,0.8,0.8,0.8,0.8,1.8,1.8,1.8,1.8),
                                    labels=c('h2', 'h21(PHE)', 'h22(PED)', 'h23(HYB)','h24(PED+HYB)', 'c*h21(PHE)', 'c*h22(PED)', 'c*h23(HYB)', 'c*h24(PED+HYB)'),
                                    name='')+
                                ylim(1.1*min(log10(d$value[d$Iteration>5])),3.5)+
                                ylab('log(h2)') +                    
                                geom_text(aes(max(d$Iteration),log10(H2),label = paste0('optim=log(',H2,')'), vjust = -1, hjust=+0.7))+
                                geom_hline(yintercept=log10(H2), linetype='dashed', color='grey25', size=1.0)+
                                geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                scale_color_manual(name="", values=c("red","green","blue"))+
                                theme(axis.text=element_text(size=14),
                                      axis.title=element_text(size=16,face="bold"),
                                      legend.title = element_text(size=16),
                                      legend.text = element_text(size=16),
                                      legend.position="bottom",
                                      legend.key.width = unit(1.9, "cm"),
                                      legend.key.height = unit(0.7, "cm"),
                                      legend.box="vertical", legend.margin=margin())+
                                guides(color     = guide_legend(order = 0),
                                       linetype  = guide_legend(order = 1),
                                       size      = guide_legend(order = 1))
        png(paste0(PATH,'eh2_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(eh2p);dev.off()
        
        #--------------------------------------------
        # Plot estim_sigHYB for all iterations
        #--------------------------------------------
        esigHYBhist<-           x_hist[,c('sigHYB', 'estim_sigHYB', 'estim_sigHYBc')]
        esigHYBhist$Iteration<- 1:nrow(esigHYBhist)
        esigHYBhist$sigHYB<-    esigHYBhist$sigHYB^2
        d <-                    melt(esigHYBhist, id.vars="Iteration")
        
        esigHYBp <-             ggplot(d, aes(Iteration, value, linetype=variable)) +
                                    geom_line(aes(size=variable))+
                                    scale_linetype_manual(values=c(1,3,5), 
                                        labels=c(expression(paste(sigma^2,'HYB')),expression(paste(sigma^2, 'HYB_estim')),expression(paste('c*',sigma^2, 'HYB_estim'))), 
                                        name='')+
                                    scale_size_manual(values=c(0.8,0.8,1.8),
                                        labels=c(expression(paste(sigma^2,'HYB')),expression(paste(sigma^2, 'HYB_estim')),expression(paste('c*',sigma^2, 'HYB_estim'))),
                                        name='')+
                                    #ylim(0,0.5)+
                                    ylim(0,1.1*max(d$value[d$Iteration>20]))+
                                    ylab(expression(paste(sigma^2,'_HYB'))) +
                                    #geom_text(aes(max(d$Iteration),0.4^2 ,label = 'optim=0.4^2', vjust = -1, hjust=+1.0))+
                                    #geom_hline(yintercept=0.4^2, linetype='dashed', color='grey25', size=1.0)+                                
                                    geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                    scale_color_manual(name="", values=c("red","green","blue"))+
                                    theme(axis.text=element_text(size=14),
                                          axis.title=element_text(size=16,face="bold"),
                                          legend.title = element_text(size=16), 
                                          legend.text = element_text(size=16),
                                          legend.position="bottom",
                                          legend.key.width = unit(1.9, "cm"),
                                          legend.key.height = unit(0.7, "cm"),
                                          legend.box="vertical", legend.margin=margin())+
                                    guides(color     = guide_legend(order = 0),
                                           linetype  = guide_legend(order = 1),
                                           size      = guide_legend(order = 1))    
        png(paste0(PATH,'esigHYB_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(esigHYBp);dev.off()
        
        #--------------------------------------------
        # Plot estim_sigSNP for all iterations
        #--------------------------------------------
        esigSNPhist<-           x_hist[,c('sigSNP', 'estim_sigSNP', 'estim_sigSNPc')]
        esigSNPhist$Iteration<- 1:nrow(esigSNPhist)
        esigSNPhist$sigSNP<-    esigSNPhist$sigSNP^2
        d <-                    melt(esigSNPhist, id.vars="Iteration")
  
        esigSNPp <-             ggplot(d, aes(Iteration, value, linetype=variable)) + geom_line(aes(size=variable))+
                                    scale_linetype_manual(values=c(1,3,5), 
                                        labels=c(expression(paste(sigma^2,'SNP')),expression(paste(sigma^2, 'SNP_estim')),expression(paste('c*',sigma^2, 'SNP_estim'))),
                                        name='')+
                                    scale_size_manual(values=c(0.8,0.8,1.8),
                                        labels=c(expression(paste(sigma^2,'SNP')),expression(paste(sigma^2, 'SNP_estim')),expression(paste('c*',sigma^2, 'SNP_estim'))),
                                        name='')+
                                    #ylim(0,0.0003)+
                                    ylim(0,1.1*max(d$value[d$Iteration>20]))+
                                    ylab(expression(paste(sigma^2,'_SNP'))) +
                                    #geom_text(aes(max(d$Iteration),0.015^2,label = 'optim=0.015^2', vjust = -1, hjust=+1.0))+
                                    #geom_hline(yintercept=0.015^2, linetype='dashed', color='grey25', size=1.0)+
                                    geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                    scale_color_manual(name="", values=c("red","green","blue"))+
                                    theme(axis.text=element_text(size=14),
                                          axis.title=element_text(size=16,face="bold"),
                                          legend.title = element_text(size=16),
                                          legend.text = element_text(size=16),
                                          legend.position="bottom",
                                          legend.key.width = unit(1.9, "cm"),
                                          legend.key.height = unit(0.7, "cm"),
                                          legend.box="vertical", legend.margin=margin())+
                                    guides(color     = guide_legend(order = 0),
                                           linetype  = guide_legend(order = 1),
                                           size      = guide_legend(order = 1))
        png(paste0(PATH,'esigSNP_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(esigSNPp);dev.off()
    } else {
        #--------------------------------------------
        # Plot r for all iterations
        #--------------------------------------------
        rhist<-                  x_hist[,c('r', 'Rep')]
        rhist$Iteration<-        1:nrow(rhist)
    
        rp <-                   ggplot(rhist, aes(Iteration, r))+
                                    geom_line(size=0.8)+
                                    #ylim(0,1)+
                                    ylim(0,1.1*max(rhist$r[rhist$Iteration>20]))+
                                    ylab(expression(paste(rho))) +
                                    geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                    scale_color_manual(name="", values=c("red","green"))+
                                    theme(axis.text=element_text(size=14),
                                          axis.title=element_text(size=16,face="bold"))
        png(paste0(PATH,'r_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(rp);dev.off()

        #--------------------------------------------
        # Plot h2 for all iterations
        #--------------------------------------------
        h2hist<-                    x_hist[,c('h2', 'Rep')]
        h2hist$Iteration<-          1:nrow(h2hist)
        h2hist$h2<-                 h2hist$h2^2
    
        h2p <-                  ggplot(h2hist, aes(Iteration, h2))+
                                    geom_line(size=0.8) +
                                    #ylim(0,2)+
                                    ylim(0,1.1*max(h2hist$h2[h2hist$Iteration>20]))+
                                    geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                    scale_color_manual(name="", values=c("red","green"))+
                                    ylab('h2') +
                                    theme(axis.text=element_text(size=14),
                                          axis.title=element_text(size=16,face="bold"))   
        png(paste0(PATH,'h2_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(h2p);dev.off()
    
        #--------------------------------------------
        # Plot sigHYB for all iterations
        #--------------------------------------------
        sigHYBhist<-            x_hist[,c('sigHYB', 'Rep')]
        sigHYBhist$Iteration<-  1:nrow(sigHYBhist)
        sigHYBhist$sigHYB<-     sigHYBhist$sigHYB^2
  
        sigHYBp <-              ggplot(sigHYBhist, aes(Iteration, sigHYB))+
                                    geom_line(size=0.8) +
                                    ylim(0,1.1*max(sigHYBhist$sigHYB[sigHYBhist$Iteration>5]))+
                                    geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                    scale_color_manual(name="", values=c("red","green"))+
                                    ylab(expression(paste(sigma^2,'_HYB'))) +
                                    theme(axis.text=element_text(size=14),
                                          axis.title=element_text(size=16,face="bold")) 
        png(paste0(PATH,'sigHYB_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(sigHYBp);dev.off()
    
        #--------------------------------------------
        # Plot sigSNP for all iterations
        #--------------------------------------------
        sigSNPhist<-            x_hist[,c('sigSNP', 'Rep')]
        sigSNPhist$Iteration<-  1:nrow(sigSNPhist)
        sigSNPhist$sigSNP<-     sigSNPhist$sigSNP^2
  
        sigSNPp <-              ggplot(sigSNPhist, aes(Iteration, sigSNP))+
                                    geom_line(size=0.8) +
                                    ylim(0,1.1*max(sigSNPhist$sigSNP[sigSNPhist$Iteration>5]))+
                                    geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                    geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                    scale_color_manual(name="", values=c("red","green"))+
                                    ylab(expression(paste(sigma^2,'_SNP'))) +
                                    theme(axis.text=element_text(size=14),
                                          axis.title=element_text(size=16,face="bold")) 
        png(paste0(PATH,'sigSNP_histPlot',SUFFIX,'.png'),width=12*300,height=7*300,res=300);plot(sigSNPp);dev.off()
    }
    }
}   