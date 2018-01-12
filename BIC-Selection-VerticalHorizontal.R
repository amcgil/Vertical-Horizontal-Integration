

VerticalHorizontalIntegration_BIC <- function(X,M,K,n,p,pcore,L,rho,lmbda1,lmbda2,lmbda3,maxiter1,maxiter2,maxiter3,tol1,tol2,tol3,tolerances,eps){
  
  N1=length(lmbda1)
  N2=length(lmbda2)
  N3=length(lmbda3)
  thetas=vector("list",N1*N2*N3)
  BIC=counts=llk=array(dim=c(N1,N2,N3))
  cnt=1
  counts.all=array(dim=c(N1,N2,N3,M,K))
  
  Xstd=X
  for(m in 1:M){
    for (k in 1:K){
      for(j in 1:p[m]){
        Xstd[[k + K*(m-1)]][,j] = X[[k + K*(m-1)]][,j]-mean(X[[k + K*(m-1)]][,j])
      }
    }
  }

  L=L+eps*diag(K)
  Lsqrt=sqrtm(L)
  Lsqrt.inv=solve(Lsqrt)
  
  temp=diag(K) + L
  temp.inv=solve(temp)
  
  thetas.zeros=list()
  for (m in 1:M){
    for (k in 1:K){
      thetas.zeros[[k + (m-1)*K]]=matrix(0,p[m],p[m])
    }
  }
  
  cnt=1
  index=array(dim=c(N1,N2,N3))
  for (i in 1:N1){
    for (j in 1:N2){
      for (l in 1:N3){
        index[i,j,l]=cnt
        cnt=cnt+1
      }
    }
  }
  
  for (i in 1:N1){
    for (j in 1:N2){
      for (l in 1:N3){
        lmbda_grid=c(lmbda1[i],lmbda2[j],lmbda3[l])
        print(lmbda_grid)
        if (i>1){
          thetas0=list()
          for (m in 1:M){
            for (k in 1:K){
              ind=index[i-1,j,l]
              thetas0[[k + (m-1)*K]]=thetas[[ind]][[k + (m-1)*K]]
            }
          }
          
          res.thetas=ADMM(Xstd, M, K, p,pcore, L, rho,lmbda_grid,maxiter1,maxiter2,maxiter3, tol1,tol2,tol3,tolerances,init=TRUE,thetas0,
               Lsqrt,Lsqrt.inv,t(Lsqrt),temp.inv)
          rm(thetas0)         
          ind=index[i,j,l]
          Rt=list()
          for (m in 1:M){
            for (k in 1:K){
              Rt[[k + (m-1)*K]]=res.thetas[[k + (m-1)*K]]
            }
          }
          thetas[[ind]]=Rt
          rm(Rt)
        } else {
          
          tolerances=c(seq(1e-4,1e-5,len=20))
          thetas.init=ADMM(Xstd, M, K, p,pcore, L, rho,lmbda_grid,maxiter1=50,maxiter2=1000,maxiter3=1000, tol1,tol2,tol3,tolerances,init=FALSE,thetas.zeros,
                          Lsqrt,Lsqrt.inv,t(Lsqrt),temp.inv)
          
          res.thetas=ADMM(Xstd, M, K, p,pcore, L, rho,lmbda_grid,maxiter1,maxiter2,maxiter3, tol1,tol2,tol3,tolerances,init=TRUE,thetas.init,
                   Lsqrt,Lsqrt.inv,t(Lsqrt),temp.inv)
          rm(thetas.init)    
          ind=index[i,j,l]
          Rt=list()
          for (m in 1:M){
            for (k in 1:K){
              Rt[[k + (m-1)*K]]=res.thetas[[k + (m-1)*K]]
            }
          }
          thetas[[ind]]=Rt
          rm(Rt)
        }
      
        BIC[i,j,l]=counts[i,j,l]=llk[i,j,l]=0
        for (m in 1:M){
          for (k in 1:K){
            Scov=cov(Xstd[[k + (m-1)*K]])	
            The=thetas[[ind]][[k + (m-1)*K]]
            The[abs(The)<1e-5]=0
            upper.entries=The[upper.tri(The)==TRUE]
            ct=length(upper.entries[abs(upper.entries)>0])
            counts.all[i,j,l,m,k]=ct
            counts[i,j,l]=counts[i,j,l]+ct
            llk[i,j,l]=llk[i,j,l]+n*determinant(The,logarithm=TRUE)$modulus[[1]] - n*tr(Scov%*%The)
            if (ct>1){
              BIC[i,j,l] = BIC[i,j,l] - n*determinant(The,logarithm=TRUE)$modulus[[1]] + n*tr(Scov%*%The) + log(n)*ct
            } else {
              BIC[i,j,l] = Inf
            }
          } # End of k loop
        } # End of m loop
      } # End of l loop
    } # End of j loop
  } # End of i loop
  
  id.BIC=which(BIC==min(BIC), arr.ind=TRUE)
  if (dim(id.BIC)[1]>1){
    id.BIC=id.BIC[1,]
  }
  
  id1=id.BIC[1]
  id2=id.BIC[2]
  id3=id.BIC[3]
  ind=index[id1,id2,id3]
  opt.tune <- c(lmbda1[id.BIC[1]],lmbda2[id.BIC[2]],lmbda3[id.BIC[3]])
  thetas.hat=thetas[[ind]]
  
  return(list(tune=opt.tune,thetas.hat=thetas.hat,BIC=BIC,counts=counts,full_counts=counts.all,thetas=thetas,BICalt=BIC.alt))
  
} # End of BIC_VH function