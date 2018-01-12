// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
cube MatrixArrayMult(NumericMatrix Mt, NumericVector Art, int pt){
  int Kt = Mt.nrow();	
  mat Mat(Mt.begin(), Kt, Kt, false);
 cube Ar(Art.begin(),pt,pt,Kt,false);
cube Lt = zeros<cube>(pt,pt,Kt);
vec v = zeros<vec>(1);
  for (int k1=0; k1<Kt; k1++){
    for (int k2=0; k2<Kt; k2++){
      mat Art = Ar.slice(k2);
      Lt.slice(k1) = Lt.slice(k1) + arma::as_scalar(Mat(k1,k2))*Art;
    }
  }
return(Lt);
}
  

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double loglik(List Thetas, List Slist, int N){
  double llk=0.0;
  
  for (int m=0; m<N; m++){
  	  mat Thl=Thetas(m);
  	  mat Sl=Slist(m);
  	  mat inn =Sl*Thl;
  double val;
double sign;

log_det(val, sign, Thl);
llk = llk + val - trace(inn);
  }
  return(llk);

}
  


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
mat soft_thres(NumericMatrix mat1, double lm) {
    int pr = mat1.nrow();
    int pc = mat1.ncol();
    mat out = zeros<mat>(pr,pc);
      for(int i = 0; i<pr; i++) {
      	for (int j=0; j<pc; j++){
        out(i,j) = std::max(fabs(mat1(i,j))-lm, 0.0);
        if (mat1(i,j)<0){
        	out(i,j)=-out(i,j);
        }
        if (mat1(i,j)==0){
        	out(i,j)=0;
        }
      }
      }
      return out;

  }


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
 mat soft_thres_square(NumericMatrix mat1, double lm) {
    int p = mat1.nrow();
    mat out =zeros<mat>(p,p);
      for(int i = 0; i<p; i++) {
      	for (int j=i; j<p; j++){
        out(i,j) = std::max(fabs(mat1(i,j))-lm, 0.0);
        if (mat1(i,j)<0){
        	out(i,j)=-out(i,j);
        }
        if (mat1(i,j)==0){
        	out(i,j)=0;
        }
        out(j,i) = out(i,j);
      }
      out(i,i)=mat1(i,i);
      }
      return out;

  }
  
  
  
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
mat is_greater(NumericMatrix Matr, int pd, double lm){
	mat ineq = zeros<mat>(pd,pd);
	for (int i=0; i<pd; i++){
		for (int j=i; j<pd; j++){
			if (Matr(i,j)>lm){
						ineq(i,j)=1;	
						ineq(j,i)=ineq(i,j);
			}
		}
	}
	return(ineq);
}
 
  
  
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List sgl_pen(List Q, int V, int pd, double rho, double lm1, double lm2){
	double lam1=lm1/rho;
	double lam2=lm2/rho;
	
	List softQ(V);
	for (int v=0; v<V; v++){
		softQ(v)=soft_thres_square(Q(v),lam1);
	}
	
	mat normsoftQ = zeros<mat>(pd,pd);
	for (int v=0; v<V; v++){
		mat Sq = softQ(v);
		normsoftQ = normsoftQ + square(Sq);
	}
	normsoftQ = sqrt(normsoftQ);
	mat notshrunk = is_greater(wrap(normsoftQ),pd,lam2);
	normsoftQ = normsoftQ + (1-notshrunk);

	List out(V);
	for (int v=0; v<V; v++){
		mat Sq = softQ(v);
		out(v)=Sq%(1-lam2/normsoftQ);
		mat Vt1=out(v);
		out(v)=Vt1%notshrunk;
		mat Vt2 = out(v);
		Vt2.diag()=Sq.diag();
		out(v)=Vt2;
	}
	return(out);
}


  

  // [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List ADMM_VH(List Xstd, int M, int K, NumericVector pv, int pcore, NumericMatrix L, double rho, NumericVector lmbda, int maxiter1,int maxiter2, int maxiter3, double tol1,double tol2,double tol3,NumericVector tolerances,bool init,List thetas_init,
NumericMatrix Lsqrt, NumericMatrix LsqrtI, NumericMatrix LsqrtT,NumericMatrix temp_inv){
	
  int G=7;
	int s=1; 
	bool conv = false;
	List S(M*K);
	List A(M*K);
	List B(M*K);
	List C(M*K);
	List C1(M*K);
	List C2(M*K);	
	List U(M*K);
	int N = tolerances.size();
	vec same_vector = zeros<vec>(G);
	int cb=1;
	mat index = zeros<mat>(M,K);
	for (int m=0; m<M; m++){
		for (int k=0; k<K; k++){
			index(m,k)=k + K*m;
		}
	}
	
	mat rA = zeros<mat>(M,K);
	mat rB = zeros<mat>(M,K);
	mat rC = zeros<mat>(M,K);
	mat rD = zeros<mat>(M,K);
	
	for (int m=0; m<M; m++){
		int pm = pv(m);
		for (int k=0; k<K; k++){
			mat Xs = Xstd(k + K*m);
			S(k + K*m)=cov(Xs);
			A(k + K*m)=zeros<mat>(pm,pm);
			B(k + K*m)=zeros<mat>(pm,pm);
			C(k + K*m)=zeros<mat>(pm,pm);
			C1(k + K*m)=zeros<mat>(pm,pm);		
			C2(k + K*m)=zeros<mat>(pm,pm);
			U(k + K*m)=zeros<mat>(pm,pm);													
		}
	}
	
	if (init==true){
        for (int m=0; m<M; m++){
            for (int k=0; k<K; k++){
                mat Ti = thetas_init(k+K*m);
                A(k+K*m)=Ti;
                B(k+K*m)=Ti;
            }
        }
	}
    
    vec vp = zeros<vec>(1);
    vec vt=zeros<vec>(1);
    
	while (s<maxiter1 && !conv){
    
        
    // FIRST ADMM
    List Bold(M*K);
        for (int m=0; m<M; m++){
            for (int k=0; k<K; k++){
                mat Bt = B(k+K*m);
                Bold(k+K*m)=Bt;
            }
        }
        
    List U1(M*K);
        for (int m=0; m<M; m++){
            int pm = pv(m);
            for (int k=0; k<K; k++){
                C(k + K*m)=zeros<mat>(pm,pm);
                U1(k + K*m)=zeros<mat>(pm,pm);
            }
        }
        
		
		int s1=1;
		bool conv1=false;
		
		while (s1<maxiter2 && !conv1){
			
            List Cold(M*K);
            for (int m=0; m<M; m++){
                for (int k=0; k<K; k++){
                    mat Ct = C(k+K*m);
                    Cold(k+K*m)=Ct;
                }
            }
			
		for (int m=0; m<M; m++){
			int pm = pv(m);
		for (int k=0; k<K; k++){
		  int id = k + K*m;
			mat Bmk = B(id); 
			mat Cmk = C(id); 
			mat Umk = U(id); 
			mat U1mk = U1(id);
			mat Smk = S(id); 
			mat Ei = rho*(Bmk + Cmk - Umk - U1mk) - Smk;
			mat Q;
			vec lmb;
			arma::eig_sym(lmb,Q, Ei);
			vec Lmb = 0.25*(lmb + sqrt(square(lmb) + 8*rho))/rho;
			mat Dl = diagmat(Lmb);
			A(k + K*m)=Q*Dl*Q.t();		
			if (m>0){
			  Bmk.resize(pm,pm);
			  Cmk.resize(pm,pm);
			  Umk.resize(pm,pm);			
			  U1mk.resize(pm,pm);
			  Smk.resize(pm,pm);		
			  Ei.resize(pm,pm);	
			  Q.resize(pm,pm);						
			}	
		}
		}
		
		for (int m=0; m<M; m++){
			int pm = pv(m);
			if (pcore<pm){
			for (int k=0; k<K; k++){
				mat Au = A(k + K*m);
				mat U1u = U1(k + K*m);
				mat temp = Au + U1u;
				mat temp1_cut = temp(span(0,pcore-1),span(pcore, pm-1));
				mat temp2_cut = temp(span(pcore,pm-1),span(0, pcore-1));
				mat temp3_cut = temp(span(pcore,pm-1),span(pcore,pm-1));
				mat v1 = soft_thres(wrap(temp1_cut),lmbda(0)/rho);
				mat v2 = soft_thres(wrap(temp2_cut),lmbda(0)/rho);
				mat v3 = soft_thres_square(wrap(temp3_cut),lmbda(0)/rho);
				mat Ctemp = zeros<mat>(pm,pm);
				Ctemp(span(0,pcore-1),span(pcore, pm-1))=v1;
				Ctemp(span(pcore,pm-1),span(0,pcore-1))=v2;
				Ctemp(span(pcore,pm-1),span(pcore, pm-1))=v3;
				C(k+m*K)=Ctemp;
				if (m>0){
				  U1u.resize(pm,pm);
				  temp.resize(pm,pm);
				  temp1_cut.resize(pcore,pm-pcore);			
				  temp2_cut.resize(pm-pcore,pcore);
				  temp3_cut.resize(pm-pcore,pm-pcore);
				  v1.resize(pcore,pm-pcore);	
				  v2.resize(pm-pcore,pcore);
				  v3.resize(pm-pcore,pm-pcore);
				  Ctemp.resize(pm,pm);
				}	
			}
			}
		}

            for (int k=0; k<K; k++){
                vec ind = index.col(k);
                List Q(M);
                for (int m=0; m<M; m++) {
                    int id = ind(m);
                    mat Ad = A(id);
                    mat U1d = U1(id);
                    mat tmp = Ad+U1d;
                    Q(m)=tmp(span(0,pcore-1),span(0,pcore-1));
                }
                List Qlist=sgl_pen(Q,M,pcore,rho,lmbda(0),lmbda(1));
                for (int m=0; m<M; m++){
                	mat Ctemp = C(k+K*m);
                    mat Atemp = A(k+K*m);
                    mat U1temp = U1(k+K*m);
                	mat Qm = Qlist(m);
                   Ctemp(span(0,pcore-1),span(0,pcore-1))=Qm;
                    Ctemp.diag()=Atemp.diag()+U1temp.diag();
                 C(k+K*m)=Ctemp;
                }
            }
		
            for (int m=0; m<M; m++){
                for (int k=0; k<K; k++){
                    mat U1mk = U1(k+K*m);
                    mat Amk = A(k+K*m);
                    mat Cmk = C(k+K*m);
                    U1(k+K*m)=U1mk+Amk-Cmk;
                }
            }
            
			double stp1 = 0.0;
			double stp2 = 0.0;
			for (int m=0; m<M; m++){
                int pm = pv(m);
				for (int k=0; k<K; k++){
					mat Amk = A(k + K*m);
					mat Cmk = C(k + K*m);
					mat Cmk_old = Cold(k + K*m);
                    mat diffAC = Amk - Cmk;
                    mat diffC = Cmk - Cmk_old;
					stp1 = stp1 + norm(diffAC,"fro");
					stp2 = stp2 + norm(diffC,"fro");
				}
			}

			double crit1 = stp1/(M*K);
			double crit2 = stp2/(M*K);
			vec crt=zeros<vec>(2);
			crt(0)=crit1;
			crt(1)=crit2;
			double crit = max(crt);
			
			if (crit<tol2){
				conv1=true;
				break;
			} else {
				s1=s1+1;
			}
	
			if (s1==maxiter2){
				conv1=true;
				break;
			}
		}
		// END OF FIRST ADMM

    
    // BEGINNING OF SECOND ADMM
    
		int s2=1;
		bool conv2 = false;
		vec critv = zeros<vec>(2);
        
        List U2(M*K);
        for (int m=0; m<M; m++){
            int pm = pv(m);
            for (int k=0; k<K; k++){
                C(k + K*m)=zeros<mat>(pm,pm);
                U2(k + K*m)=zeros<mat>(pm,pm);
            }
        }
    
		
		while (s2<maxiter3 && !conv2){
		  
		  double lam3=lmbda(2)/rho;
		  for (int m=0; m<M; m++){
		    int pm=pv(m);
		    rowvec indr = index.row(m);
		    cube Cr = zeros<cube>(pm,pm,K);
		    cube U2r = zeros<cube>(pm,pm,K);
		    if (m>0){
		      Cr.resize(pm,pm,K);
		      U2r.resize(pm,pm,K);
		    }
		    for (int k=0; k<K;k++){
		      int id=indr(k);
		      mat Ct = C(id);
		      mat U2t = U2(id);
		      if (m>0){
		        Ct.resize(pm,pm);
		        U2t.resize(pm,pm);
		      }
		      Cr.slice(k)=Ct;
		      U2r.slice(k)=U2t;
		    }
            
		    cube Lc=MatrixArrayMult(Lsqrt,wrap(Cr),pm);
		    cube Lu=MatrixArrayMult(LsqrtI,wrap(U2r),pm);
              
		    cube t1=zeros<cube>(pm,pm,K);
		    cube t2=zeros<cube>(pm,pm,K);			
		    mat normtemp = zeros<mat>(pm,pm);
		    for (int k=0; k<K; k++){
		      t1.slice(k)=Lc.slice(k) - U2r.slice(k);
		      t2.slice(k)=Cr.slice(k) - Lu.slice(k);
		      mat tls = t1.slice(k);
		      normtemp = normtemp+ square(tls);
		    }
		    normtemp= sqrt(normtemp);

		    mat notshrunk = is_greater(wrap(normtemp),pm,lam3);
		    normtemp = normtemp + (1-notshrunk);
              
		    List out(K);
		    for (int k=0; k<K; k++){
		      mat t2m = t2.slice(k);
		      out(k)=t2m%(1-lam3/normtemp);
		      mat Vt1=out(k);
		      out(k)=Vt1%notshrunk;
		      mat Vt2 = out(k);
		      Vt2.diag()=t2m.diag();
		      out(k)=Vt2;
		      int id=indr(k);
		      B(id)=out(k);
		    }
              
		  }
		  
            List Cold(M*K);
		          for (int m=0; m<M; m++){
                      for (int k=0; k<K; k++){
                          mat Ct = C(k+K*m);
                          Cold(k+K*m)=Ct;
                      }
                  }
            
		  for (int m=0; m<M; m++){
		    int pm=pv(m);
		    rowvec indr = index.row(m);
		    cube Ar = zeros<cube>(pm,pm,K);
		    cube Br = zeros<cube>(pm,pm,K);
		    cube Ur = zeros<cube>(pm,pm,K);
		    cube U2r = zeros<cube>(pm,pm,K);
              
		    if (m>0){
		      Ar.resize(pm,pm,K);
		      Br.resize(pm,pm,K);
		      Ur.resize(pm,pm,K);
		      U2r.resize(pm,pm,K);
		    }
		    
		    for (int k=0; k<K; k++){
		      int id=indr(k);
		      mat At = A(id);
		      mat Bt = B(id);
		      mat Ut = U(id);
		      mat U2t = U2(id);
		      
		      if (m>0){
		        At.resize(pm,pm);
		        Bt.resize(pm,pm);
		        Ut.resize(pm,pm);
		        U2t.resize(pm,pm);
		      }
		      
		      Ar.slice(k)=At;
		      Br.slice(k)=Bt;
		      Ur.slice(k)=Ut;
		      U2r.slice(k)=U2t;
		    }
		    
		    cube Lb=MatrixArrayMult(L,wrap(Br),pm);
              
		    cube Lu2=MatrixArrayMult(LsqrtT,wrap(U2r),pm);
              
		    cube SumArray=Ar+Ur+Lb+Lu2;
              
		    cube Cr=MatrixArrayMult(temp_inv,wrap(SumArray),pm);
		    
		    for (int k=0; k<K; k++){
		      int id = indr(k);
		      C(id)=Cr.slice(k);
		    }
              
		  }
            

		  for (int m=0; m<M; m++){
		    int pm=pv(m);
		    rowvec indr = index.row(m);
		    cube Br = zeros<cube>(pm,pm,K);
		    cube Cr = zeros<cube>(pm,pm,K);
		    cube U2r = zeros<cube>(pm,pm,K);	
		    
		    if (m>0){
		      Br.resize(pm,pm,K);
		      Cr.resize(pm,pm,K);
		      U2r.resize(pm,pm,K);
		    }
		    
		    for (int k=0; k<K; k++){
		      int id=indr(k);
		      mat Bt = B(id);
		      mat Ct = C(id);
		      mat U2t = U2(id);
		      
		      if (m>0){
		        Bt.resize(pm,pm);
		        Ct.resize(pm,pm);
		        U2t.resize(pm,pm);
		      }
		      
		      Br.slice(k)=Bt;
		      Cr.slice(k)=Ct;
		      U2r.slice(k)=U2t;
		    }
		    cube DiffAr = zeros<cube>(pm,pm,K);
      
		    
		    for (int k=0; k<K; k++){
		      DiffAr.slice(k)=Br.slice(k)-Cr.slice(k);
            }
            
		    U2r = U2r + MatrixArrayMult(Lsqrt,wrap(DiffAr),pm);
		    
		    for (int k=0; k<K; k++){
		      int id = indr(k);
		      U2(id)=U2r.slice(k);
		    }						
		  }
    
		  
		  mat rB = zeros<mat>(M,K);
		  mat rC = zeros<mat>(M,K);
		  for (int m=0; m<M; m++){
              int pm = pv(m);
		    for (int k=0; k<K; k++){
		      mat Bd = B(k + m*K);
		      mat Cd = C(k + m*K);
		      mat Cd_old = Cold(k + m*K);
		      mat diff = Bd - Cd;
		      mat gr = Cd - Cd_old;
		      rB(m,k)=norm(diff,"fro");
		      rC(m,k)=norm(gr,"fro");
		    }
		  }
            
		  
		  critv(0)=arma::as_scalar(accu(rB));
		  critv(1)=arma::as_scalar(accu(rC));
		  double Mx=max(critv)/(M*K);

		  if (Mx<tol3){
		    conv2=true;
		    break;
		  } else {
		    s2=s2+1;
		    conv2=false;
		  }	
		  
		  if (s2==maxiter3){
		    conv2=true;
		    break;
		  }	
		  
      
            
		}
    // END of SECOND ADMM
        
        for (int m=0; m<M; m++){
            for (int k=0; k<K; k++){
                mat Umk = U(k+K*m);
                mat Amk = A(k+K*m);
                mat Bmk = B(k+K*m);
                U(k+K*m)=Umk+Amk-Bmk;
            }
        }
     
	double stp1 = 0.0;
	double stp2 = 0.0;
		for (int m=0; m<M; m++){
            int pm = pv(m);
			for (int k=0; k<K; k++){
				mat Amk = A(k + K*m);
				mat Bmk = B(k + K*m);
				mat Bmk_old = Bold(k + K*m);
                mat diffAB = Amk - Bmk;
                mat diffB = Bmk - Bmk_old;
				stp1 = stp1 + norm(diffAB,"fro");
				stp2 = stp2 + norm(diffB,"fro");
			}
		}
	double crit1 = stp1/(M*K);
	double crit2 = stp2/(M*K);
	vec crt=zeros<vec>(2);
	crt(0)=crit1;
	crt(1)=crit2;
	double crit = max(crt);	

        rowvec vc = zeros<rowvec>(2);
        vc(0)=s;
        vc(1)=crit;
        vc.print();
        
        	int nedges = 0;
		int nedges_old = 0;

        if (crit<0.05){
           bool same = false;
            for (int m=0; m<M; m++){
                for (int k=0; k<K; k++){
                    mat Bmk = B(k+K*m);
                    mat Bmk_old = Bold(k+K*m);
                    same = approx_equal(Bmk,Bmk_old,"absdiff",1e-3);
                    if (same==false){
                        break;
                    }
                }
            }


 		for (int m=0; m<M; m++){
                for (int k=0; k<K; k++){
                    mat Bmk = B(k+K*m);
                    mat Bmk_old = Bold(k+K*m);
			  int pm = pv(m);
     			for (int i=0; i<pm; i++){
			for (int j=i; j<pm; j++){
				if (fabs(Bmk(i,j))<1e-4){
				nedges = nedges+1;
			}
				if (fabs(Bmk_old(i,j))<1e-4){
				nedges_old = nedges_old+1;
			}
			}
			}

                }
            }

            
            if (same==true && cb<=G){
                same_vector(cb-1)=1;
                cb=cb+1;
            }
            if (cb>G){
                cb=1;
            }
            
            rowvec vc2 = zeros<rowvec>(G+2);
            for (int i=0; i<G; i++){
                vc2(i)=same_vector(i);
            }
			vc2(G)=nedges;
			vc2(G+1)=nedges_old;
            vc2.print();
            

        }
      
        /*mat Aprint = A(0);
         mat Aprint_cut = Aprint(span(0,3),span(0,3));
         for (int i=0; i<4; i++){
         for (int j=0; j<4; j++){
         if (fabs(Aprint_cut(i,j))<1e-3){
         Aprint_cut(i,j)=0;
         Aprint_cut(j,i)=Aprint_cut(i,j);
         }
         }
         }
         Aprint_cut.print();*/
    
	if (s<=N){
	tol2=tolerances(s-1);
	tol3=tol2;
	}
        
        if (crit<1e-3 && arma::as_scalar(accu(same_vector))==G && nedges==nedges_old){
            conv=true;
            break;
        }
	
	if (crit<tol1){
		conv=true;
		break;
	} else {
		s=s+1;
	}
		
		
	if (s==maxiter1){
		conv=true;
		break;
	}
	}

	return(A);
		
}

    
  
  
  
  
  
  
  
  
