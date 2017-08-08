CEFF<- function(winn){
  # returns the 0.99 cef for each window
  # win<-c("ID",.25,.05) #choose "TR", "ID", "PR", or "QS" 
  # with any parameters that apply
  
  if(winn[1]=="TR")
  {cc<-as.numeric(winn[2])
  cef=cc}
  
  if(winn[1]=="QS")
  {
    b<-as.numeric(winn[2])
    cc<-as.numeric(winn[3])
    cef=cc}
  
  if(winn[1]=="PR")
  {
    cc<-as.numeric(win[2])
    cef=cc}
  
  if(winn[1]=="ID")
  {
    b<-as.numeric(winn[2])
    cc<-as.numeric(winn[3])
    xx=c(1:1000)/1000
    yy=xx
    for (i in 1:1000){
      yy[i]=LambdaIDbc(xx[i],b,cc)}
    cef=max(xx[yy>.99])}
  
  cef}





LambdaIDbc <- function(x,b=1/4,c)
{
  # A function to compute Infinitely Differentible LambdaTRbc
  # x<-input of the function
  # b<-shape parameter, b > 0
  # c<-c in (0,1] --good values are in [.05, .15]
  
  y<-abs(x)
  
  
  if ( y   <= c)
  {
    val<-1
  } 
  else if ( (c < y) && (y < 1) )
  {
    denom<-(y-1)^2
    numer<- -b*exp(-b/ ((y-c)^2) )
    z<-numer/denom
    val<-exp(z)
  }
  else if (y >= 1)
  {
    val<-0
  }
  
  
  return(val)
  
}








###################################################
CmpGamma <- function(Vt, short=T)
{
  # A function to compute Autocovariance matrix sequence Gamma
  # Vt<-input d-variate time series
  # len<-T<-Length of Time Series 
  # short=T only computes the autocovariance up to lag  2*sqrt(len)
  dimen<-dim(Vt)
  d<-dimen[1]
  len<-dimen[2]
  LEN=len-1
  if (short==T) {LEN=floor(2*sqrt(len)) }
  
  Gamma<-array(c(1:(d*d*len))*0,dim=c(d,d,len))
  
  for(j in 1:LEN)
  {
    temp <- matrix(data=0, nrow=d, ncol=d) 
    for (k in 1: (len-j))            
    {
      temp<-temp + Vt[,k] %*% t(Vt[,(k+j)])                
    }
    temp<-temp/len
    Gamma[,,j]<-temp
  }
  
  # special case for j<-0, it is placed at Gamma(len) i.e Gamma(T)
  
  j<-0
  temp <- matrix(data=0, nrow=d, ncol=d) 
  for (k in 1: (len-j))            
  {
    temp<-temp + Vt[,k] %*% t(Vt[,(k+j)])                
  }
  temp<-temp/len
  Gamma[,,len]<-temp
  
  
  return (Gamma)
  
}

EvalGammaJ <- function(Gamma,j)
{
  # A function to evaluate Gamma(j) 
  # for any value of j from negative through zero to positive integer
  # Gamma<-Computed Gamma series for j<-0 to T-1 its dimension is (d,d,T)
  # j<-index of Gamma
  # len<-T<-Length of Time Series 
  
  
  dimen<-dim(Gamma)
  len<-dimen[3]
  
  if (j == 0)
  {
    GJ <- (Gamma[,,len])
  }
  else if (j < 0)
  {
    k<--j
    temp <- Gamma[,,k]
    GJ <- (t(temp))
  }
  else
    GJ <- (Gamma[,,j])
  
  return (GJ) 
  
}

CmpRho <- function(Gamma, short=T)
{
  # short=T only computes the autocovariance up to lag 2*sqrt(len)
  
  dimen<-dim(Gamma)
  d<-dimen[1]
  len<-dimen[3]   
  
  LEN=len-1
  if (short==T) {LEN=floor(2*sqrt(len)) }
  
  Rho<-array(c(1:(d*d*len))*0,dim=c(d,d,len))
  
  
  Gam0 <- EvalGammaJ(Gamma,0)
  norm<-sqrt(diag(Gam0) %*%t(diag(Gam0)))
  
  for (m in 1: LEN)   
  { Rho[,,m]<-Gamma[,,m]/norm   
  }
  Rho[,,len]<-Gamma[,,len]/norm 
  
  return (Rho)
  
  
}

EvalRhoM <- function(Rho,m)
{
  # A function to evaluate Rho(m) 
  # for any value of nonnegative integer m 
  # Rho<-Computed correlogram/cross-correlogram for m<-0 to T-1 its dimension is (d,d,T)
  # m<-index of Rho
  # len<-T<-Length of Time Series 
  
  dimen<-dim(Rho)
  len<-dimen[3]
  
  if (m == 0)
  {
    RJ <- (Rho[,,len])
  }
  else
    RJ <- (Rho[,,m])
  
  return (RJ) 
  
  
}



############ 

CmpHatS <- function(Rho,short=T)
{
  # A function to compute Bandwidth Matrix HatS from Rho
  # Rho<-Estimated correlogram/cross-correlogram Matrix of dimension (d,d,T)
  # d<-Dimension of each input vector V(t)
  # len<-T<-Length of Time Series 
  # cef<-cef value given
  
  dimen<-dim(Rho)
  d<-dimen[1]
  len<-dimen[3]
  
  
  
  HatS<-matrix(nrow=d, ncol=d)
  HatQ <-matrix(nrow=d, ncol=d)
  
  C0<-2
  
  # KT<-max(5,sqrt(log10(T)))  , here KT is an integer valued function of T
  
  KT<-max(5,sqrt(log(len,10)))
  
  # val<-C0*sqrt((log10(T))/T)  , val is the value to compare with
  
  LEN=len-KT-1
  if (short==T) {LEN=min(floor(2*sqrt(len)),len-KT-1) }
  
  
  p<-log10(len)
  p<-p/len
  p<-sqrt(p)
  val<-C0*p
  
  
  for(j in 1:d)
  {
    for(k in 1:d) 
    { 
      for (q in 1: LEN)
      {
        found<-1
        for (m in 1: KT)
        {
          RhoJKQM <- abs(EvalRhoM(Rho,(q+m)))
          if ( RhoJKQM[j,k] >= val) 
          {
            found<-0
            break 
          }
        } # for m         
        
        if ( found == 1)
        {
          HatQ[j,k]<-ifelse(is.na(q),floor(sqrt(len)),q)
          break
        }    # if (found == 1)     
      }  # for q
      
      # Add this line to avoid NA values
      HatQ[j,k]<-ifelse(is.na(HatQ[j,k]),floor(sqrt(len)),HatQ[j,k])
      
    } # for k
  } # for j
  
  
  
  for(j in 1:d)
  {
    for(k in 1:d) 
    { 
      if (j == k)
      {
        x<-HatQ[j,k]/cef
        x<-ceiling(x)
        x<-max(x,1)
        HatS[j,k]<-x 
      } else
      {
        hatq<-max(HatQ[j,k],HatQ[k,j])
        x<-hatq/cef
        x<-ceiling(x)
        x<-max(x,1)
        HatS[j,k]<-x
        HatS[k,j]<-x
        
      }
    } # for k
  } # for j
  
  
  
  return (HatS)
  
}
##############







