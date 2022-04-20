################################################################################
#################################meOPERATORS.r##################################
################################################################################
# R-Script/function to define operators used in meBLUP MML
#
#
# Himmelbauer, 14.11.2021
##################################################################################
library(data.table)
library(Matrix)
rnd <- function(x) trunc(x+sign(x)*0.5)

##########################################
# Input section
##########################################
    nSNP<<-  config$nSNP
    prt<<-   config$prt
    na<<-    config$na   #number of animals
    nap<<-   config$nap  #number of animals with phenotypes
    nf<<-    config$nf   #number of fixed effects 
    nas<<-   config$nas  #numberof animals with genotypes
    ns<<-    config$ns   #number of SNPs
    m<<-     config$m
    n<<-     config$n
    q<-      config$q
    nup<<-   as.numeric(classDef[,5]==2)+as.numeric(classDef[,6]==2)
   
##########################################
# Define Operators
##########################################
rnd <- function(x) trunc(x+sign(x)*0.5)
#-------------------------------------------------
#Covariance matrix
#-------------------------------------------------
    # Factors pi
    pip<-        function(x) rep(1/(x[(nf+na+nSNP+2)]^2)-1, nap)     #pi for phenotype equations (used in C)
    pipd<-       function(x) rep(-1, nap)                            #pi for phenotype equations (used in C')
    pigh<-       function(x) ((nup[A$ac]+2)/4)*(x[(nf+na+nSNP+3)]^2) #pi for hybrid    equations (used in C and C')
    pigr<-       function(x) ((nup[A$ac]+2)/4)*(x[(nf+na+nSNP+2)]^2) #pi for hybrid    equations (used in C and C')
    pis<-        function(x) rep(1, nSNP)                            #pi for snp       equations (used in C and C')
    
    sdlcC<-function(x){ #sum(diag(chol(Cx))) - !!! only for univariate Version!!!
        h2<-      x[(nf+na+nSNP+2)]^2
        sigHYB2<- x[(nf+na+nSNP+3)]^2
        sigSNP<-  x[(nf+na+nSNP+4)]^2
        sig2PHE<- (1-h2)
        
        nup<-    as.numeric(classDef[,5]==2)+as.numeric(classDef[,6]==2)
        cD<<-    sum(log(sqrt(sig2PHE*rep(1, nap))))
        cG<<-    sum(log(sqrt(((nup[A$ac]+2)/4))*sqrt(h2)*sqrt(sigHYB2)))
        cM<<-    sum(log(sqrt(sigSNP) *rep(1, nSNP)))
       
        sdlcCx<-  sum(cD, cG, cM)
        return(sdlcCx) 
    }

    Cix<-function(x){ # whole C'
        h2<-      x[(nf+na+nSNP+2)]^2
        sigHYB2<- x[(nf+na+nSNP+3)]^2
        sigSNP<-  x[(nf+na+nSNP+4)]^2
        sig2PHE<- (1-h2)
        
        nup<-    as.numeric(classDef[,5]==2)+as.numeric(classDef[,6]==2)
        cD<<-    sig2PHE*rep(1, nap)
        cG<<-    ((nup[A$ac]+2)/4)*h2*sigHYB2
        cM<<-    sigSNP *rep(1, nSNP)
        
        Ci<-     Diagonal((nap+na+nSNP),c(1/cD, 1/cG,1/cM))
        return(Ci)
    }

    Cp<-function(x){ #C for phenotype equations
        h2<-         x[(nf+na+nSNP+2)]^2
        sig2PHE<-    (1-h2)
        #Cphe<-      diag(sig2PHE*rep(1, nap), nap, nap)
        Cphe<-       Diagonal(nap, sig2PHE*rep(1, nap))
        return(Cphe)
    }

    Cg<-function(x){ #C for hybrid equations
        h2<-         x[(nf+na+nSNP+2)]^2
        sigHYB2<-    x[(nf+na+nSNP+3)]^2
        nup<-        as.numeric(classDef[,5]==2)+as.numeric(classDef[,6]==2)
        #Cgen<-      diag(((nup[A$ac]+2)/4)*h2, na, na)
        Cgen<-       Diagonal(na,((nup[A$ac]+2)/4)*h2*sigHYB2)
        return(Cgen)
    }

    Cs<-function(x){ #C for snp equations
        sigSNP<-     x[(nf+na+nSNP+4)]
        #Csnp<-      diag(sigSNP^2*rep(1, nSNP), nSNP, nSNP)
        Csnp<-       Diagonal(nSNP,sigSNP^2*rep(1, nSNP))
        return(Csnp)
    }    
    
#-------------------------------------------------
#Functions used in f_MML and gradient
#-------------------------------------------------
    X <-         function(u) A$X%*%u 
    Z<-          function(u) u[A$ap]
    P<-          function(u) (A$PA%*%u)
    zGEN<<-      rep(0,na)
    ZS<-         function(v) { zs<-zGEN; zs[A$as]<-S%*%v-1*sum(v) ; return(zs) } #Aenderung
    
    Ax<-function(x) { #operator multiplying by A
        beta<-   x[1:nf]
        u<-      x[(nf+1):(nf+na)]
        v<-      x[(nf+na+1):(nf+na+nSNP)]
        r<-      x[(nf+na+nSNP+1)]
        if (!all(is.finite(x))) {
            stop("Bad input for Ax!")
        }
        Rd<<-    rep(1, na)
        if (nSNP>0) { Rd[A$as]<<-   r }
        # equations
        # D=PHE  - nap  phenotype equations:     zD=X*beta+Z*u   
        # G=HYB  - na   hybrid equations:        zG=Z_S*v-u+P*u
        # M=REGM - nSNP regularizing equations:  zM=Iv
        if (nSNP>0) {ZSv<- ZS(v)} else {ZSv<-rep(0, na)}
        if (nSNP>0) {
        z<-      rbind((X(beta)+Z(u)), (ZSv-u+Rd*P(u)), t(t(v)))
        } else {
        z<-      rbind((X(beta)+Z(u)), (ZSv-u+Rd*P(u)))
        }
        return(z)
    }
    
    r<-          function(x) Ax(x)-Fy           #whole r
    s<-          function(x) Cix(x)%*%r(x)      #whole s
     
    Ep<-         function(x) { #E for phenotype equations
                    sCp<-     solve(Cp(x))
                    spx<-     sCp%*%rp(x)
                    spx2<-    spx^2 #tcrossprod(spx)
                    ep<-      diag((2*q)*sCp)-spx2
                    return(list(ep, spx))
                    }
    Eg<-         function(x) { #E for hybrid equations
                    sCg<-     solve(Cg(x))
                    sgx<-     sCg%*%rg(x)
                    sgx2<-    sgx^2 #tcrossprod(sgx) 
                    eg<-      diag((2*q)*sCg)-sgx2
                    return(list(eg, sgx))
                    }
    Es<-         function(x) {  #E for snp equations
                    sCs<-     solve(Cs(x))
                    ssx<-     sCs%*%rs(x)
                    ssx2<-    ssx^2 #tcrossprod(ssx)
                    es<-      diag((2*q)*sCs)-ssx2
                    return(list(es, ssx))
                    }
            
    Lgamma<-     function(x,g) { #L_gamma
        if (g==1) { L<-    x[(nf+na+nSNP+2)] }
        if (g==2) { L<-    x[(nf+na+nSNP+3)] }
        if (g==3) { L<-    x[(nf+na+nSNP+4)] }
        if (g>3)  { stop("Bad input for gamma!") }
        return(L)
    }
 
#-------------------------------------------------
#Function to calcualte EBVs for test animals
#-------------------------------------------------
    candsBV<-    function(x){
        beta<-   x[1:nf]
        u<-      x[(nf+1):(nf+na)]
        v<-      x[(nf+na+1):(nf+na+nSNP)]
        r<-      x[(nf+na+nSNP+1)]
        if (!all(is.finite(x))) {
            stop("Bad input for Ax!")
        }
        cands<-  A$can
        cangt<-  A$can%in%A$ascan
        
        Rdcan<-         data.table(ID=A$can, resid=rep(1, length(A$can)))
        Rdcan$resid[Rdcan$ID%in%cangt]<-r
        Rdcan<-         Rdcan$resid
        # hybrid equations:        u=S*v+r*PA*u
        uGEN<-        rep(0,length(A$can))
        if (nSNP>0) {
        uGEN[cangt]<- A$Scan%*%v-1*sum(v)
        }
        uPA<-         Rdcan*(A$PAcan%*%u)
        uCan<-        data.table(ID=A$can, EBV=as.vector(uGEN+uPA))
        return(uCan)    
    }

#-------------------------------------------------
#Function to derive scaling vector for starting values
#-------------------------------------------------
    meSCALE=function(x) {
        b<-           x[(1)     :(nf)]
        u<-           x[(nf+1)   :(nf+na)]
        v<-           x[(nf+na+1):(nf+na+nSNP)]
        r<-           x[(nf+na+nSNP+1)]
        h<-           x[(nf+na+nSNP+2)]
        sigHYB<-      x[(nf+na+nSNP+3)]
        sigSNP<-      x[(nf+na+nSNP+4)]
    
    
        sb<-          max(b)
        su<-          max(u)
        sv<-          max(v)
        sr<-          max(r)
        sh<-          max(h)
        ssr<-         max(sigHYB)
        sss<-         max(sigSNP)
    
        if (nSNP>0) {
        scl<-         c(rep(sb,nf), rep(su, na), rep(sv, nSNP),sr, sh, ssr, sss)
        } else {
        scl<-         c(rep(sb,nf), rep(su, na), rep(sv, nSNP),sr, sh, ssr)
        }
        return(scl) 
    }



