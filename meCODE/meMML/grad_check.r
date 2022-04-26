########################################################################
##############################gradcheck.r###############################
########################################################################
# function [relacc,bad]=gradcheck(getfg,x,prt);
# gradient check for evaluation routine 
#    results depend on random choices 
#    repeat with several x for higher reliability
#
# getfg          # function handle
#                # [f,g]=getfg(x) is supposed to compute 
#                # for a column vector of size n 
#                # a function value f=f(x) and its gradient g=f'(x)^T
# x              # point of evaluation 
#                #   (should be generic, not very close to stationary)
# prt            # print level (1 = print something)
#
# relacc         # estimated relative accuracy of gradient
# bad            # list of bad gradient components
#



gradcheck<-function(getf, getg, x, prt, idx){

eps<-    2.220446049250313e-16

ni<-     15 # number of trial evaluations
qi<-     rep(0, ni)

n<-      length(x)
bad<-    rep(0, n)
x[x=0]<- 0.001

# function value and gradient at starting point
f0<-    getf(x)
df<-    getg(x)

# random direction with sensible scale
p<-     (1+runif(n,0,1))*abs(x)
p[p=0]  <-abs(x[p=0]) # enforce nonzero components

#Use idx to define which gradients are checked and in which order
if(missing(idx)) { idx<- order(abs(getg(x))) } 
kidx<- c(0, idx)


for (k in kidx) {
  if (k>0) { 
    # check k-th component of gradient separately
    if (prt) {
        print(paste0('check ',k,'th component of gradient separately'))
    }
    p<-   rep(0, n)
    p[k]<-abs(x[k])
  }
  gTp<-   t(df)%*%t(t(p))
  for (i in 1:ni) {
    alp<-  eps*10^i
    fi<-   getf(x+alp*p)
    qi[i]<-(fi-f0)/(alp*gTp)
    if (is.nan(qi[i])) {
    qi[i]<-0
    warning("nan")
    #print(paste0(i,k,'nan'))
    }    
  }
  acc<-    abs(qi-1)
  # most of the first few components should be <<1
  if (prt & k==0) {
    print('most of the first few numbers in the following list')
    print('should be <<1:')
    print(acc)
    print(' ')
    
  }
  good<-   0
  maxgood<-0
  for (i in 1:ni) {
    #print(acc)
    if (acc[i]<=0.1) {
        good<-  good+1
    } else {
        good<-  0
    }
    if (good>maxgood) {
      maxgood<-   good
      igood<-     i
    }
  }
  if (maxgood>2) {
    # directional derivative seems ok
    if (k==0) {
      ind<-   igood+1-maxgood:igood
      relacc<-min(acc[ind])
      break
    }
    } else {
    # directional derivative seems erroneous
    if (k==0) {
      relacc=min(acc[acc>0.1])
    } else {
      bad[k]<- 1
    }
  }

}
# find list of bad components
bad<-  which(bad==1)
return(bad)
}
