################################################################################
#################################meGRAD.r#######################################
################################################################################
# R-Script/function for gradient of f_mml
#
#-------------------------------------------
#USAGE:
# grad(x)
#   x....point at which the gradient should be calculated
#
#--------------------------------------------
#
# Himmelbauer, 14.11.2021 | last update 19.04.2022
##################################################################################

library(data.table)
library(Matrix)
rnd <- function(x) trunc(x+sign(x)*0.5)

#Gradient
grad<-  function(x){
    verbose=1
    x<-             x*scl
    
    u<-             x[(nf+1):(nf+na)]
    v<-             x[(nf+na+1):(nf+na+nSNP)]
    rd<-            x[(nf+na+nSNP+1)]
    Rd<-            rep(1, na)
    Rd[A$as]<-      rd
       
    #fixe Effekte (db)
    db<-            list(A$X, Matrix(nrow=na, ncol=nf, data=0, sparse=TRUE), Matrix(nrow=nSNP, ncol=nf, data=0, sparse=TRUE))
        
    #animal effects (du)
    Zap<-           Matrix(nrow=nap, ncol=na, data=0, sparse=TRUE); diag(Zap[1:nap, A$ap])<-1
    du<-            list(Zap, as((Diagonal(na, (-1))+A$PA%*%Diagonal(length(Rd), Rd)), "sparseMatrix"), Matrix(nrow=nSNP, ncol=na, data=0, sparse=TRUE))    
        
    #SNP Effekte (dv)
    ZSNP<-          matrix(0, na, nSNP)
    ZSNP[A$as,]<-   S-1
    dv<-            list(Matrix(nrow=nap, ncol=nSNP, data=0, sparse=TRUE), ZSNP, Diagonal(nSNP))

    #R (dR)
    Pu2<-           as(P(u),"sparseMatrix")
    Pu<-            Matrix(nrow=na, ncol=1, data=0, sparse=TRUE)
    Pu[A$as,]<-     Pu2[A$as,]      
    dr<-            list(Matrix(nrow=nap, ncol=1, data=0, sparse=TRUE), Pu, Matrix(nrow=nSNP, ncol=1, data=0, sparse=TRUE))
    
    #Egamma (1-3)
    tmp<-           Ep(x); Epx<-tmp[[1]]; spx<-tmp[[2]]
    tmp<-           Eg(x); Egx<-tmp[[1]]; sgx<-tmp[[2]]
    tmp<-           Es(x); Esx<-tmp[[1]]; ssx<-tmp[[2]]
    Eg1<-           sum(pipd(x)*Epx)+sum(pigh(x)*Egx)
    Eg2<-           sum(pigr(x)*Egx)
    Eg3<-           sum(pis(x) *Esx)
    
    #s
    if (nSNP>0) {
    st<-            list(t(spx), t(sgx), t(ssx))
    } else {
    st<-            list(t(spx), t(sgx))
    }
    
    u0<-            as.matrix(t(rep(0,na)))
    v0<-            as.matrix(t(rep(0,nSNP)))
    
    # Gradient
    if (nSNP>0) {
        gr<-as.vector(cbind(
            st[[1]]%*%db[[1]]+st[[2]]%*%db[[2]]+st[[3]]%*%db[[3]],        #beta   f'MML=s^T*my'       
            st[[1]]%*%du[[1]]+st[[2]]%*%du[[2]]+st[[3]]%*%du[[3]],        #u      f'MML=s^T*my'  
            st[[1]]%*%dv[[1]]+st[[2]]%*%dv[[2]]+st[[3]]%*%dv[[3]],        #v      f'MML=s^T*my'
            st[[1]]%*%dr[[1]]+st[[2]]%*%dr[[2]]+st[[3]]%*%dr[[3]],        #rho    f'MML=s^T*my'
            Eg1%*% Lgamma(x,1),                                           #h2      
            Eg2%*% Lgamma(x,2),                                           #sigHYB       
            Eg3%*% Lgamma(x,3)                                            #sigSNP    
    ))
    } else {
        gr<-as.vector(cbind(
            st[[1]]%*%db[[1]]+st[[2]]%*%db[[2]]+st[[3]]%*%db[[3]],        #beta   f'MML=s^T*my'       
            st[[1]]%*%du[[1]]+st[[2]]%*%du[[2]]+st[[3]]%*%du[[3]],        #u      f'MML=s^T*my'  
            st[[1]]%*%dr[[1]]+st[[2]]%*%dr[[2]]+st[[3]]%*%dr[[3]],        #rho    f'MML=s^T*my'
            Eg1%*% Lgamma(x,1),                                           #h2      
            Eg2%*% Lgamma(x,2),                                           #sigHYB          
    ))
    }
    
    gr<-            gr*scl
    
    if (verbose) {
        fwrite(data.table(it=gradcount, gnorm=max(abs(gr)), idx=which.max(abs(gr))), 'GRAD.out', col.names=F, row.names=F,quote=F,sep=' ',dec='.', append=TRUE)
        gradcount<<-gradcount+1
    }
    
    return(gr)
}
    
#Gradient conv. EBV  
grad_EBV<-          function(x){
    verbose=1
    x<-             x*scl
    
    u<-             x[(nf+1):(nf+na)]
    v<-             x[(nf+na+1):(nf+na+nSNP)]
    rd<-            x[(nf+na+nSNP+1)]
    Rd<-            rep(1, na)
    Rd[A$as]<-      rd
       
    #fixe Effekte (db)
    db<-             list(A$X, Matrix(nrow=na, ncol=nf, data=0, sparse=TRUE), Matrix(nrow=nSNP, ncol=nf, data=0, sparse=TRUE))
        
    #animal effects (du)
    Zap<-            Matrix(nrow=nap, ncol=na, data=0, sparse=TRUE); diag(Zap[1:nap, A$ap])<-1
    du<-             list(Zap, as((Diagonal(na, (-1))+A$PA%*%Diagonal(length(Rd), Rd)), "sparseMatrix"), Matrix(nrow=nSNP, ncol=na, data=0, sparse=TRUE))    
        
    #SNP Effekte (dv)
    ZSNP<-          matrix(0, na, nSNP)
    ZSNP[A$as,]<-   S-1
    dv<-            list(Matrix(nrow=nap, ncol=nSNP, data=0, sparse=TRUE), ZSNP, Diagonal(nSNP))

    #R (dR)
    Pu2<-            as(P(u),"sparseMatrix")
    Pu<-             Matrix(nrow=na, ncol=1, data=0, sparse=TRUE)
    dr<-             list(Matrix(nrow=nap, ncol=1, data=0, sparse=TRUE), Pu, Matrix(nrow=nSNP, ncol=1, data=0, sparse=TRUE))
    
    #Egamma (1-3)
    tmp<-     Ep(x); Epx<-tmp[[1]]; spx<-tmp[[2]]
    tmp<-     Eg(x); Egx<-tmp[[1]]; sgx<-tmp[[2]]
    
    Eg1<-     sum(pipd(x)*Epx)+sum(pigh(x)*Egx)
    Eg2<-     sum(pigr(x)*Egx)
    
    #s
    st<-            list(t(spx), t(sgx))
    
    u0<-            as.matrix(t(rep(0,na)))
    v0<-            as.matrix(t(rep(0,nSNP)))
    
    #f'MML=s^T*my'
    gr<-as.vector(cbind(
        st[[1]]%*%db[[1]]+st[[2]]%*%db[[2]],                    #beta          
        st[[1]]%*%du[[1]]+st[[2]]%*%du[[2]],                    #u       
        0,                                                      #r
        0,                                                      #h2      
        0                                                       #sigHYB           
    ))
    
    gr<-   gr*scl
    
    if (verbose) {
        fwrite(data.table(it=gradcount, gnorm=max(abs(gr)), idx=which.max(abs(gr))), 'GRAD.out', col.names=F, row.names=F,quote=F,sep=' ',dec='.', append=TRUE)
        gradcount<<-gradcount+1
    }
    
    return(gr)
}
    

#Gradient GEBV1
grad_GEBV1<-          function(x){
    verbose=1
    x<-               x*scl
    
    u<-               x[(nf+1):(nf+na)]
    v<-               x[(nf+na+1):(nf+na+nSNP)]
    rd<-              x[(nf+na+nSNP+1)]
    Rd<-              rep(1, na)
    Rd[A$as]<-        rd
       
    #fixe Effekte (db)
    db<-              list(A$X, Matrix(nrow=na, ncol=nf, data=0, sparse=TRUE), Matrix(nrow=nSNP, ncol=nf, data=0, sparse=TRUE))
        
    #animal effects (du)
    Zap<-             Matrix(nrow=nap, ncol=na, data=0, sparse=TRUE); diag(Zap[1:nap, A$ap])<-1
    du<-              list(Zap, as((Diagonal(na, (-1))+A$PA%*%Diagonal(length(Rd), Rd)), "sparseMatrix"), Matrix(nrow=nSNP, ncol=na, data=0, sparse=TRUE))    
        
    #SNP Effekte (dv)
    ZSNP<-           matrix(0, na, nSNP)
    ZSNP[A$as,]<-    S-1
    dv<-             list(Matrix(nrow=nap, ncol=nSNP, data=0, sparse=TRUE), ZSNP, Diagonal(nSNP))

    #R (dR)
    Pu2<-            as(P(u),"sparseMatrix")
    Pu<-             Matrix(nrow=na, ncol=1, data=0, sparse=TRUE)
    Pu[A$as,]<-      Pu2[A$as,] 
    dr<-             list(Matrix(nrow=nap, ncol=1, data=0, sparse=TRUE), Pu, Matrix(nrow=nSNP, ncol=1, data=0, sparse=TRUE))
    
    #Egamma (1-3)
    tmp<-            Ep(x); Epx<-tmp[[1]]; spx<-tmp[[2]]
    tmp<-            Eg(x); Egx<-tmp[[1]]; sgx<-tmp[[2]]
    tmp<-            Es(x); Esx<-tmp[[1]]; ssx<-tmp[[2]]
    
    Eg1<-            sum(pipd(x)*Epx)+sum(pigh(x)*Egx)
    Eg2<-            sum(pigr(x)*Egx)
    Eg3<-            sum(pis(x) *Esx)

    #s
    st<-             list(t(spx), t(sgx), t(ssx))
    
    u0<-             as.matrix(t(rep(0,na)))
    v0<-             as.matrix(t(rep(0,nSNP)))
    
    #f'MML=s^T*my'
    gr<-as.vector(cbind(
        0,                                                             #beta          
        u0,                                                            #u       
        st[[1]]%*%dv[[1]]+st[[2]]%*%dv[[2]]+st[[3]]%*%dv[[3]],         #v  
        0,                                                             #r
        0,                                                             #h2      
        0,                                                             #sigHYB       
        Eg3%*% Lgamma(x,3)                                             #sigSNP    
    ))
    
    gr<-             gr*scl
    
    if (verbose) {
        fwrite(data.table(it=gradcount, gnorm=max(abs(gr)), idx=which.max(abs(gr))), 'GRAD.out', col.names=F, row.names=F,quote=F,sep=' ',dec='.', append=TRUE)
        gradcount<<-gradcount+1
    }
    
    return(gr)
}


#Gradient GEBV2
grad_GEBV2<-         function(x){
    verbose=1
    x<-              x*scl
    
    u<-              x[(nf+1):(nf+na)]
    v<-              x[(nf+na+1):(nf+na+nSNP)]
    rd<-             x[(nf+na+nSNP+1)]
    Rd<-             rep(1, na)
    Rd[A$as]<-       rd
       
    #fixe Effekte (db)
    db<-             list(A$X, Matrix(nrow=na, ncol=nf, data=0, sparse=TRUE), Matrix(nrow=nSNP, ncol=nf, data=0, sparse=TRUE))
        
    #animal effects (du)
    Zap<-            Matrix(nrow=nap, ncol=na, data=0, sparse=TRUE); diag(Zap[1:nap, A$ap])<-1
    du<-             list(Zap, as((Diagonal(na, (-1))+A$PA%*%Diagonal(length(Rd), Rd)), "sparseMatrix"), Matrix(nrow=nSNP, ncol=na, data=0, sparse=TRUE))    
        
    #SNP Effekte (dv)
    ZSNP<-           matrix(0, na, nSNP)
    ZSNP[A$as,]<-    S-1
    dv<-             list(Matrix(nrow=nap, ncol=nSNP, data=0, sparse=TRUE), ZSNP, Diagonal(nSNP))

    #R (dR)
    Pu2<-            as(P(u),"sparseMatrix")
    Pu<-             Matrix(nrow=na, ncol=1, data=0, sparse=TRUE)
    if (nSNP>0) {
    Pu[A$as,]<-      Pu2[A$as,] 
    }
    dr<-             list(Matrix(nrow=nap, ncol=1, data=0, sparse=TRUE), Pu, Matrix(nrow=nSNP, ncol=1, data=0, sparse=TRUE))
    
    #Egamma (1-3)
    tmp<-            Ep(x); Epx<-tmp[[1]]; spx<-tmp[[2]]
    tmp<-            Eg(x); Egx<-tmp[[1]]; sgx<-tmp[[2]]
    tmp<-            Es(x); Esx<-tmp[[1]]; ssx<-tmp[[2]]
    
    Eg1<-            sum(pipd(x)*Epx)+sum(pigh(x)*Egx)
    Eg2<-            sum(pigr(x)*Egx)
    Eg3<-            sum(pis(x) *Esx)

    
    #s
    st<-             list(t(spx), t(sgx), t(ssx))

    u0<-             as.matrix(t(rep(0,na)))
    v0<-             as.matrix(t(rep(0,nSNP)))
    
    #f'MML=s^T*my'
    gr<-             as.vector(cbind(
                        st[[1]]%*%db[[1]]+st[[2]]%*%db[[2]]+st[[3]]%*%db[[3]],         #beta          
                        st[[1]]%*%du[[1]]+st[[2]]%*%du[[2]]+st[[3]]%*%du[[3]],         #u       
                        st[[1]]%*%dv[[1]]+st[[2]]%*%dv[[2]]+st[[3]]%*%dv[[3]],         #v  
                        0,                                                             #r
                        0,                                                             #h2      
                        0,                                                             #sigHYB       
                        0                                                              #sigSNP    
    ))
   
    gr<-             gr*scl
    
    if (verbose) {
        fwrite(data.table(it=gradcount, gnorm=max(abs(gr)), idx=which.max(abs(gr))), 'GRAD.out', col.names=F, row.names=F,quote=F,sep=' ',dec='.', append=TRUE)
        gradcount<<- gradcount+1
    }
    
    return(gr)
}












