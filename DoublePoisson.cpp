#include <RcppArmadillo.h>
#include <cmath>




using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector my_fun_d(NumericVector y, NumericVector g){		// works perfectly fine
   
   int n = y.size();
   NumericVector res(n);
   for(int i=0; i<n; ++i){
		if (y[i] == 0) {
			res[i] = 2*g[i];
		} else {
			res[i] = 2*(y[i]*log(y[i])-y[i]-y[i]*g[i]+exp(g[i]));
		}
	}
   
    return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double my_fun_di(double y, double g){		
   
	double res;
	if (y == 0) {
		res = 2*g;
	} else {
		res = 2*(y*log(y)-y-y*g+exp(g));
	}
	
    return res;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List my_fun_PoisM(NumericVector xseq, NumericVector x, NumericVector y, List startv, const double kappa, NumericVector gammai){
	
	int n = y.size();
	int nt = xseq.size();
		
	double pi = 3.14159265358979323846264338327950;
	double beskappa = Rf_bessel_i(kappa,0,1);
	
	List startvv = Rcpp::clone(startv);
	NumericVector beta0 = startvv[0];
	NumericVector beta1 = startvv[1];
	
	NumericMatrix x_t(n,nt);
	NumericMatrix s(n,nt);
	NumericMatrix ets(n,nt);
	NumericMatrix kkxtM(n,nt);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			x_t(i,j) = x[i]-xseq[j];
			s(i,j) = sin(x_t(i,j));
			ets(i,j) = exp(beta0[j]+beta1[j]*s(i,j));
			kkxtM(i,j) = (exp(kappa * cos (x_t(i,j)))/(beskappa*2*pi))/gammai[i];
		}
	}
	
	
	
	NumericVector bn0(nt);
	NumericVector bn1(nt);
	NumericVector bn2(nt);
	for(int j=0; j<nt; ++j){
		
		bn0[j] = sum(kkxtM(_,j)*ets(_,j));
		bn1[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j));
		bn2[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j)*s(_,j));

	}
	
	
	NumericMatrix bjti0(n,nt);
	NumericMatrix bjti1(n,nt);
	NumericMatrix K(n,nt);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			bjti0(i,j) = kkxtM(i,j)*(bn2[j]-s(i,j)*bn1[j]);
			bjti1(i,j) = kkxtM(i,j)*(s(i,j)*bn0[j]-bn1[j]);	
			K(i,j) = y[i]-ets(i,j);			
		}
	}
	

	for(int j=0; j<nt; ++j){
		beta0[j] = sum(bjti0(_,j)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta0[j];
		beta1[j] = sum(bjti1(_,j)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta1[j];
	}
	
	List L = List::create(beta0, beta1);
	return L;
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List my_fun_PoisG(NumericVector xseq, NumericVector x, NumericVector y, List startv2, const double nu, NumericVector gi){
	
	int n = y.size();
	int nt = xseq.size();
		
	double pi = 3.14159265358979323846264338327950;
	double besnu = Rf_bessel_i(nu,0,1);

	List startvv2 = Rcpp::clone(startv2);	
	NumericVector beta0G = startvv2[0];
	NumericVector beta1G = startvv2[1];
	
	NumericMatrix x_t(n,nt);
	NumericMatrix s(n,nt);
	NumericMatrix tsG(n,nt);
	NumericMatrix etsG(n,nt);
	NumericMatrix kkxtG(n,nt);
	NumericVector di(n);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			x_t(i,j) = x[i]-xseq[j];
			s(i,j) = sin(x_t(i,j));
			tsG(i,j) = beta0G[j]+beta1G[j]*s(i,j);
			etsG(i,j) = exp(tsG(i,j));                                  
			kkxtG(i,j) = exp(nu * cos (x_t(i,j)))/(besnu*2*pi);
		}
		di[i] =  my_fun_di(y[i],gi[i]);
	}
	
	NumericMatrix q(n,nt);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			q(i,j) = tsG(i,j) + di[i]/etsG(i,j) - 1;
		}
	}


	NumericVector sn0G(nt);
	NumericVector sn1G(nt);
	NumericVector sn2G(nt);
	for(int j=0; j<nt; ++j){
		sn0G[j] = sum(kkxtG(_,j));
		sn1G[j] = sum(kkxtG(_,j)*s(_,j));
		sn2G[j] = sum(kkxtG(_,j)*s(_,j)*s(_,j));
	}
	

	NumericMatrix bjti0G(nt,n);
	NumericMatrix bjti1G(nt,n);
	for(int j=0; j<nt; ++j){
		for(int i=0; i<n; ++i){
			bjti0G(j,i) = kkxtG(i,j)*(sn2G[j]-s(i,j)*sn1G[j]);
			bjti1G(j,i) = kkxtG(i,j)*(s(i,j)*sn0G[j]-sn1G[j]);		
		}
	}
	

	for(int j=0; j<nt; ++j){
		beta0G[j] = sum(bjti0G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
		beta1G[j] = sum(bjti1G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
	}
	
	List L = List::create(beta0G, beta1G);
	
	return L;
}












// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List my_fun_DoublePois(NumericVector xseq, NumericVector x, NumericVector y, NumericVector startv, NumericVector startv2, const double kappa, const double nu){
	
	int n = y.size();
	int nt = xseq.size();
	
	NumericVector psii1 = rep(startv2[0],n);
	NumericVector psii2 = rep(startv2[1],n);
	List psii = List::create(psii1,psii2);

	List starting = List::create(rep(startv[0],n),rep(startv[1],n));
	
	double pi = 3.14159265358979323846264338327950;
	double beskappa = Rf_bessel_i(kappa,0,1);
	double besnu = Rf_bessel_i(nu,0,1);
	
	
	NumericVector beta0 = rep(startv[0],nt);
	NumericVector beta1 = rep(startv[1],nt);
	
	NumericMatrix x_t(n,nt);
	NumericMatrix s(n,nt);
	NumericMatrix t_s(n,nt);
	NumericMatrix ets(n,nt);
	NumericMatrix kxtM(n,nt);
	NumericMatrix kkxtM(n,nt);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			x_t(i,j) = x[i]-xseq[j];
			s(i,j) = sin(x_t(i,j));
			t_s(i,j) = beta0[j]+beta1[j]*s(i,j);
			ets(i,j) = exp(t_s(i,j));
			kxtM(i,j) = exp(kappa * cos (x_t(i,j)))/(beskappa*2*pi);
			kkxtM(i,j) = kxtM(i,j)/exp(psii1[i]);
		}
	}
	
	NumericVector bn0(nt);
	NumericVector bn1(nt);
	NumericVector bn2(nt);
	for(int j=0; j<nt; ++j){	
		bn0[j] = sum(kkxtM(_,j)*ets(_,j));
		bn1[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j));
		bn2[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j)*s(_,j));
	}
	
	
	NumericMatrix bjti0(nt,n);
	NumericMatrix bjti1(nt,n);
	NumericMatrix K(n,nt);
	for(int j=0; j<nt; ++j){
		for(int i=0; i<n; ++i){
			bjti0(j,i) = kkxtM(i,j)*(bn2[j]-s(i,j)*bn1[j]);
			bjti1(j,i) = kkxtM(i,j)*(s(i,j)*bn0[j]-bn1[j]);		
			K(i,j) = y[i]-ets(i,j);
		}
	}
	
	for(int j=0; j<nt; ++j){
		beta0[j] = sum(bjti0(j,_)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta0[j];
		beta1[j] = sum(bjti1(j,_)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta1[j];
	}

	// Dispersion part
	NumericVector beta0G = rep(startv2[0],nt);
	NumericVector beta1G = rep(startv2[1],nt);
	
	NumericMatrix tsG(n,nt);
	NumericMatrix etsG(n,nt);
	NumericMatrix kkxtG(n,nt);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
	
			tsG(i,j) = beta0G[j]+beta1G[j]*s(i,j);
			etsG(i,j) = exp(tsG(i,j));
			kkxtG(i,j) = exp(nu * cos (x_t(i,j)))/(besnu*2*pi);
		}
	}
	
	
	List gi = my_fun_PoisM(x,x,y,starting,kappa,exp(psii1));
	NumericVector di = my_fun_d(y,gi[0]);
	

	NumericMatrix q(n,nt);
	for(int j=0; j<nt; ++j){
		q(_,j) = tsG(_,j) + di/etsG(_,j) - 1;
	}

	NumericVector sn0G(nt);
	NumericVector sn1G(nt);
	NumericVector sn2G(nt);
	for(int j=0; j<nt; ++j){
		sn0G[j] = sum(kkxtG(_,j));
		sn1G[j] = sum(kkxtG(_,j)*s(_,j));
		sn2G[j] = sum(kkxtG(_,j)*s(_,j)*s(_,j));
	}
	

	NumericMatrix bjti0G(nt,n);
	NumericMatrix bjti1G(nt,n);
	for(int j=0; j<nt; ++j){
		for(int i=0; i<n; ++i){
			bjti0G(j,i) = kkxtG(i,j)*(sn2G[j]-s(i,j)*sn1G[j]);
			bjti1G(j,i) = kkxtG(i,j)*(s(i,j)*sn0G[j]-sn1G[j]);		
		}
	}
	
	
	// beta0G<- colSums(t(bjti0G) * q)/(sn2G*sn0G-sn1G^2)
	// beta1G <- colSums(t(bjti1G) * q)/(sn2G*sn0G-sn1G^2)
	for(int j=0; j<nt; ++j){
		beta0G[j] = sum(bjti0G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
		beta1G[j] = sum(bjti1G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
	}
	
	

	
	NumericVector beta0_old;
	NumericVector beta1_old;
	NumericVector beta0G_old;
	NumericVector beta1G_old;
	List psii_old;
	List gi_old;
	double res1;
	double res2;
	double res3;
	double res4;
	NumericVector resi;
	

	
	// while loop
	double tol = 0.00001;
	double res = 50000;
	int maxit = 100;
	int numit = 0;
	while ((res>=1+tol || res <= 1-tol) & (numit < maxit)){	
		beta0_old = Rcpp::clone(beta0);
		beta1_old = Rcpp::clone(beta1);
		beta0G_old = Rcpp::clone(beta0G);
		beta1G_old = Rcpp::clone(beta1G);
		psii_old = Rcpp::clone(psii);
		
		psii = my_fun_PoisG(x,x,y,psii_old,nu,gi[0]);
		psii1 = psii[0];
		for(int i=0; i<n; ++i){
			for(int j=0; j<nt; ++j){
				kkxtM(i,j) = kxtM(i,j)/exp(psii1[i]);
				t_s(i,j) = beta0[j]+beta1[j]*s(i,j);
				ets(i,j) = exp(t_s(i,j));
			}
		}
		
		for(int j=0; j<nt; ++j){	
			bn0[j] = sum(kkxtM(_,j)*ets(_,j));
			bn1[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j));
			bn2[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j)*s(_,j));
		}
		
		for(int j=0; j<nt; ++j){
			for(int i=0; i<n; ++i){
				bjti0(j,i) = kkxtM(i,j)*(bn2[j]-s(i,j)*bn1[j]);
				bjti1(j,i) = kkxtM(i,j)*(s(i,j)*bn0[j]-bn1[j]);		
				K(i,j) = y[i]-ets(i,j);
			}
		}
		
		for(int j=0; j<nt; ++j){
			beta0[j] = sum(bjti0(j,_)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta0[j];
			beta1[j] = sum(bjti1(j,_)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta1[j];
		}


		// Dispersion part
		for(int i=0; i<n; ++i){
			for(int j=0; j<nt; ++j){
				tsG(i,j) = beta0G[j]+beta1G[j]*s(i,j);
				etsG(i,j) = exp(tsG(i,j));
			}
		}
		gi_old = gi;
		gi = my_fun_PoisM(x,x,y,gi_old,kappa,exp(psii1));
		di = my_fun_d(y,gi[0]);

		for(int j=0; j<nt; ++j){
			q(_,j) = tsG(_,j) + di/etsG(_,j) - 1;
		}
	
	
		for(int j=0; j<nt; ++j){
			sn0G[j] = sum(kkxtG(_,j));
			sn1G[j] = sum(kkxtG(_,j)*s(_,j));
			sn2G[j] = sum(kkxtG(_,j)*s(_,j)*s(_,j));
		}
		
		for(int j=0; j<nt; ++j){
			for(int i=0; i<n; ++i){
				bjti0G(j,i) = kkxtG(i,j)*(sn2G[j]-s(i,j)*sn1G[j]);
				bjti1G(j,i) = kkxtG(i,j)*(s(i,j)*sn0G[j]-sn1G[j]);		
			}
		}
		
		
		// beta0G<- colSums(t(bjti0G) * q)/(sn2G*sn0G-sn1G^2)
		// beta1G <- colSums(t(bjti1G) * q)/(sn2G*sn0G-sn1G^2)
		for(int j=0; j<nt; ++j){
			beta0G[j] = sum(bjti0G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
			beta1G[j] = sum(bjti1G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
		}
		//res<-max(c(abs(beta0G/beta0G_old),abs(beta1G/beta1G_old),abs(beta0/beta0_old),abs(beta1/beta1_old)))
		res1 = max(abs(beta0G/beta0G_old));
		res2 = max(abs(beta1G/beta1G_old));
		res3 = max(abs(beta0/beta0_old));
		res4 = max(abs(beta1/beta1_old));
		
		resi = {res1,res2,res3,res4};
		res = max(resi);

		numit = numit+1;
	}

	List L = List::create(beta0,beta0G);

	return L;
	
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List my_fun_DoublePois_j(NumericVector xseq, NumericVector x, NumericVector y, NumericVector startv, NumericVector startv2, const double kappa, const double nu, int k){
	
	int kk = k-1;
	int n = y.size();
	int nt = xseq.size();
	
	NumericVector psii1 = rep(startv2[0],n);
	NumericVector psii2 = rep(startv2[1],n);
	List psii = List::create(psii1,psii2);

	List starting = List::create(rep(startv[0],n),rep(startv[1],n));
	
	double pi = 3.14159265358979323846264338327950;
	double beskappa = Rf_bessel_i(kappa,0,1);
	double besnu = Rf_bessel_i(nu,0,1);
	
	
	NumericVector beta0 = rep(startv[0],nt);
	NumericVector beta1 = rep(startv[1],nt);
	
	NumericMatrix x_t(n,nt);
	NumericMatrix s(n,nt);
	NumericMatrix ets(n,nt);
	NumericMatrix kxtM(n,nt);
	NumericMatrix kkxtM(n,nt);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			x_t(i,j) = x[i]-xseq[j];
			s(i,j) = sin(x_t(i,j));
			ets(i,j) = exp(beta0[j]+beta1[j]*s(i,j));
			kxtM(i,j) = exp(kappa * cos (x_t(i,j)))/(beskappa*2*pi);
			kkxtM(i,j) = kxtM(i,j)/exp(psii1[i]);
		}
	}
	
	NumericVector bn0(nt);
	NumericVector bn1(nt);
	NumericVector bn2(nt);
	for(int j=0; j<nt; ++j){	
		bn0[j] = sum(kkxtM(_,j)*ets(_,j));
		bn1[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j));
		bn2[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j)*s(_,j));
	}
	
	
	NumericMatrix bjti0(n,nt);
	NumericMatrix bjti1(n,nt);
	NumericMatrix K(n,nt);
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			bjti0(i,j) = kkxtM(i,j)*(bn2[j]-s(i,j)*bn1[j]);
			bjti1(i,j) = kkxtM(i,j)*(s(i,j)*bn0[j]-bn1[j]);	
			K(i,j) = y[i]-ets(i,j);			
		}
	}
	

	for(int j=0; j<nt; ++j){
		beta0[j] = sum(bjti0(_,j)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta0[j];
		beta1[j] = sum(bjti1(_,j)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta1[j];
	}

	// Dispersion part
	NumericVector beta0G = rep(startv2[0],nt);
	NumericVector beta1G = rep(startv2[1],nt);
	
	NumericVector x_kk = Rcpp::clone(x);
	x_kk.erase(kk);
	
	NumericMatrix x_t_kk(n-1,nt);
	NumericMatrix s_kk(n-1,nt);
	NumericMatrix ts_kk(n-1,nt);
	NumericMatrix ets_kk(n-1,nt);
	
	for(int i=0; i<(n-1); ++i){
		for(int j=0; j<nt; ++j){
				x_t_kk(i,j) = x_kk[i]-xseq[j];
				s_kk(i,j) = sin(x_t_kk(i,j));
				ts_kk(i,j) = beta0[j]+beta1[j]*s_kk(i,j);
				ets_kk(i,j) = exp(ts_kk(i,j));	
		}
	}
	
	
	NumericMatrix tsG(n-1,nt);
	NumericMatrix etsG(n-1,nt);
	NumericMatrix kkxtG(n-1,nt);
	for(int i=0; i<(n-1); ++i){
		for(int j=0; j<nt; ++j){
			tsG(i,j) = beta0G[j]+beta1G[j]*s_kk(i,j);
			etsG(i,j) = exp(tsG(i,j));
			kkxtG(i,j) = exp(nu * cos (x_t_kk(i,j)))/(besnu*2*pi);
		}
	}
	
	
	List gi = my_fun_PoisM(x,x,y,starting,kappa,exp(psii1));
	NumericVector gi0 = gi[0];
	NumericVector y_kk = Rcpp::clone(y);
	y_kk.erase(kk);
	NumericVector gi0_kk = Rcpp::clone(gi0);	
	gi0_kk.erase(kk);
	NumericVector di = my_fun_d(y_kk,gi0_kk);
	

	NumericMatrix q(n-1,nt);
	for(int i=0; i<(n-1); ++i){
		for(int j=0; j<nt; ++j){
			q(i,j) = tsG(i,j) + di[i]/etsG(i,j) - 1;
		}
	}
	

	NumericVector sn0G(nt);
	NumericVector sn1G(nt);
	NumericVector sn2G(nt);
	for(int j=0; j<nt; ++j){
		sn0G[j] = sum(kkxtG(_,j));
		sn1G[j] = sum(kkxtG(_,j)*s_kk(_,j));
		sn2G[j] = sum(kkxtG(_,j)*s_kk(_,j)*s_kk(_,j));
	}
	

	NumericMatrix bjti0G(nt,n-1);
	NumericMatrix bjti1G(nt,n-1);
	for(int j=0; j<nt; ++j){
		for(int i=0; i<(n-1); ++i){
			bjti0G(j,i) = kkxtG(i,j)*(sn2G[j]-s_kk(i,j)*sn1G[j]);
			bjti1G(j,i) = kkxtG(i,j)*(s_kk(i,j)*sn0G[j]-sn1G[j]);		
		}
	}
	
	
	for(int j=0; j<nt; ++j){
		beta0G[j] = sum(bjti0G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
		beta1G[j] = sum(bjti1G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
	}
	
	
	NumericVector beta0_old;
	NumericVector beta1_old;
	NumericVector beta0G_old;
	NumericVector beta1G_old;
	List psii_old;
	List gi_old;
	double res1;
	double res2;
	double res3;
	double res4;
	NumericVector resi;
	

	
	// while loop
	double tol = 0.00001;
	double res = 50000;
	while (res>=1+tol || res <= 1-tol){
	
		beta0_old = Rcpp::clone(beta0);				
		beta1_old = Rcpp::clone(beta1);				
		beta0G_old = Rcpp::clone(beta0G);			
		beta1G_old = Rcpp::clone(beta1G);			
		psii_old = Rcpp::clone(psii);				
		
		psii = my_fun_PoisG(x,x_kk,y_kk,psii_old,nu,gi0_kk);
		psii1 = psii[0];
		
		for(int i=0; i<n; ++i){
			for(int j=0; j<nt; ++j){
				kkxtM(i,j) = kxtM(i,j)/exp(psii1[i]);
				ets(i,j) = exp(beta0[j]+beta1[j]*s(i,j));
			}
		}
		
		for(int j=0; j<nt; ++j){	
			bn0[j] = sum(kkxtM(_,j)*ets(_,j));
			bn1[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j));
			bn2[j] = sum(kkxtM(_,j)*ets(_,j)*s(_,j)*s(_,j));
		}
		
	
		for(int i=0; i<n; ++i){
			for(int j=0; j<nt; ++j){
				bjti0(i,j) = kkxtM(i,j)*(bn2[j]-s(i,j)*bn1[j]);
				bjti1(i,j) = kkxtM(i,j)*(s(i,j)*bn0[j]-bn1[j]);		
				K(i,j) = y[i]-ets(i,j);
			}
		}
		
		for(int j=0; j<nt; ++j){
			beta0[j] = sum(bjti0(_,j)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta0[j];
			beta1[j] = sum(bjti1(_,j)*K(_,j))/(bn2[j]*bn0[j]-bn1[j]*bn1[j]) + beta1[j];
		}


		// Dispersion part
		for(int i=0; i<(n-1); ++i){
			for(int j=0; j<nt; ++j){
				tsG(i,j) = beta0G[j]+beta1G[j]*s_kk(i,j);
				etsG(i,j) = exp(tsG(i,j));
			}
		}
		gi_old = gi;
		gi = my_fun_PoisM(x,x,y,gi_old,kappa,exp(psii1));
		gi0_kk = gi[0];
		gi0_kk.erase(kk);
		di = my_fun_d(y_kk,gi0_kk);

		NumericMatrix q(n-1,nt);
		for(int i=0; i<(n-1); ++i){
			for(int j=0; j<nt; ++j){
				q(i,j) = tsG(i,j) + di[i]/etsG(i,j) - 1;
			}
		}
	
		for(int j=0; j<nt; ++j){
			sn0G[j] = sum(kkxtG(_,j));
			sn1G[j] = sum(kkxtG(_,j)*s_kk(_,j));
			sn2G[j] = sum(kkxtG(_,j)*s_kk(_,j)*s_kk(_,j));
		}
		
		for(int j=0; j<nt; ++j){
			for(int i=0; i<(n-1); ++i){
				bjti0G(j,i) = kkxtG(i,j)*(sn2G[j]-s_kk(i,j)*sn1G[j]);
				bjti1G(j,i) = kkxtG(i,j)*(s_kk(i,j)*sn0G[j]-sn1G[j]);		
			}
		}
		
		for(int j=0; j<nt; ++j){
			beta0G[j] = sum(bjti0G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
			beta1G[j] = sum(bjti1G(j,_)*q(_,j))/(sn2G[j]*sn0G[j]-sn1G[j]*sn1G[j]);
		}
		
		//res<-max(c(abs(beta0G/beta0G_old),abs(beta1G/beta1G_old),abs(beta0/beta0_old),abs(beta1/beta1_old)))
		res1 = max(abs(beta0G/beta0G_old));
		res2 = max(abs(beta1G/beta1G_old));
		res3 = max(abs(beta0/beta0_old));
		res4 = max(abs(beta1/beta1_old));
		
		resi = {res1,res2,res3,res4};
		res = max(resi);

	}

	List L = List::create(beta0,beta0G);

	return L;
	
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List my_fun_DoublePois_j_v2(double xseq, NumericVector x, NumericVector y, NumericVector startv, NumericVector startv2, const double kappa, const double nu, int k){
	
	int kk = k-1;
	int n = y.size();

	
	NumericVector psii1 = rep(startv2[0],n);
	NumericVector psii2 = rep(startv2[1],n);
	List psii = List::create(psii1,psii2);

	List starting = List::create(rep(startv[0],n),rep(startv[1],n));
	
	double pi = 3.14159265358979323846264338327950;
	double beskappa = Rf_bessel_i(kappa,0,1);
	double besnu = Rf_bessel_i(nu,0,1);
	
	double beta0 = startv[0];
	double beta1 = startv[1];
	
	NumericVector x_t(n);
	NumericVector s(n);
	NumericVector t_s(n);
	NumericVector ets(n);
	NumericVector kxtM(n);
	NumericVector kkxtM(n);
	for(int i=0; i<n; ++i){
			x_t(i) = x[i]-xseq;
			s(i) = sin(x_t(i));
			t_s(i) = beta0+beta1*s(i);
			ets(i) = exp(t_s(i));
			kxtM(i) = exp(kappa * cos (x_t(i)))/(beskappa*2*pi);
			kkxtM(i) = kxtM(i)/exp(psii1[i]);
	}
	
	double bn0 = 0;
	double bn1 = 0;
	double bn2 = 0;	
	for(int i=0; i<n; ++i){
		
		bn0 += kkxtM[i]*ets[i];
		bn1 += kkxtM[i]*ets[i]*s[i];
		bn2 += kkxtM[i]*ets[i]*s[i]*s[i];
	}


	
	NumericVector bjti0(n);
	NumericVector bjti1(n);
	NumericVector K(n);

	for(int i=0; i<n; ++i){
		bjti0(i) = kkxtM(i)*(bn2-s(i)*bn1);
		bjti1(i) = kkxtM(i)*(s(i)*bn0-bn1);		
		K(i) = y[i]-ets(i);
	}

	double aux1 = 0;
	double aux2 = 0;
	for(int i=0; i<n; ++i){
		aux1 += bjti0[i]*K[i];
		aux2 += bjti1[i]*K[i];
	}

	beta0 = aux1/(bn2*bn0-bn1*bn1) + beta0;
	beta1 = aux2/(bn2*bn0-bn1*bn1) + beta1;

	// Dispersion part
	double beta0G = startv2[0];
	double beta1G = startv2[1];
	
	NumericVector x_kk = Rcpp::clone(x);
	x_kk.erase(kk);
	
	NumericVector x_t_kk(n-1);
	NumericVector s_kk(n-1);
	NumericVector ts_kk(n-1);
	NumericVector ets_kk(n-1);
	
	for(int i=0; i<(n-1); ++i){
		x_t_kk(i) = x_kk[i]-xseq;
		s_kk(i) = sin(x_t_kk(i));
		ts_kk(i) = beta0+beta1*s_kk(i);
		ets_kk(i) = exp(ts_kk(i));	
	}
	
	NumericVector tsG(n-1);
	NumericVector etsG(n-1);
	NumericVector kkxtG(n-1);
	for(int i=0; i<(n-1); ++i){
		tsG(i) = beta0G+beta1G*s_kk(i);
		etsG(i) = exp(tsG(i));
		kkxtG(i) = exp(nu * cos (x_t_kk(i)))/(besnu*2*pi);
	}
	
	List gi = my_fun_PoisM(x,x,y,starting,kappa,exp(psii1));
	NumericVector gi0 = gi[0];
	NumericVector y_kk = Rcpp::clone(y);
	y_kk.erase(kk);
	NumericVector gi0_kk = Rcpp::clone(gi0);	
	gi0_kk.erase(kk);
	NumericVector di = my_fun_d(y_kk,gi0_kk);
	
	
	NumericVector q(n-1);
	for(int i=0; i<(n-1); ++i){
		q[i] = tsG[i] + di[i]/etsG[i] - 1;
	}
	
	double sn0G = 0;
	double sn1G = 0;
	double sn2G = 0;
	for(int i=0; i<(n-1); ++i){
		sn0G += kkxtG[i];
		sn1G += kkxtG[i]*s_kk[i];
		sn2G += kkxtG[i]*s_kk[i]*s_kk[i];
	}
	
	NumericVector bjti0G(n-1);
	NumericVector bjti1G(n-1);
	for(int i=0; i<(n-1); ++i){
		bjti0G(i) = kkxtG(i)*(sn2G-s_kk(i)*sn1G);
		bjti1G(i) = kkxtG(i)*(s_kk(i)*sn0G-sn1G);		
	}

	aux1 = 0;
	aux2 = 0;
	for(int i=0; i<(n-1); ++i){
		aux1 += bjti0G[i]*q[i];
		aux2 += bjti1G[i]*q[i];	
	}
	beta0G = aux1/(sn2G*sn0G-sn1G*sn1G);
	beta1G = aux2/(sn2G*sn0G-sn1G*sn1G);

	
	
	double beta0_old;
	double beta1_old;
	double beta0G_old;
	double beta1G_old;
	List psii_old;
	List gi_old;
	NumericVector resi;
	
	// while loop
	double tol = 0.00001;
	double res = 50000;
	int maxit = 100;
	int numit = 0;
	while ((res>=1+tol || res <= 1-tol) & (numit < maxit)){
		
		beta0_old = beta0;				
		beta1_old = beta1;				
		beta0G_old = beta0G;		
		beta1G_old = beta1G;			
		psii_old = psii;			
		
		psii = my_fun_PoisG(x,x_kk,y_kk,psii_old,nu,gi0_kk);
		psii1 = psii[0];
		
		for(int i=0; i<n; ++i){
			kkxtM(i) = kxtM(i)/exp(psii1[i]);
			t_s(i) = beta0+beta1*s(i);
			ets(i) = exp(t_s(i));
		}
		
		bn0 = 0;
		bn1 = 0;
		bn2 = 0;	
		for(int i=0; i<n; ++i){
			bn0 += kkxtM[i]*ets[i];
			bn1 += kkxtM[i]*ets[i]*s[i];
			bn2 += kkxtM[i]*ets[i]*s[i]*s[i];
		}


		for(int i=0; i<n; ++i){
			bjti0(i) = kkxtM(i)*(bn2-s(i)*bn1);
			bjti1(i) = kkxtM(i)*(s(i)*bn0-bn1);		
			K(i) = y[i]-ets(i);
		}

		aux1 = 0;
		aux2 = 0;
		for(int i=0; i<n; ++i){
			aux1 += bjti0[i]*K[i];
			aux2 += bjti1[i]*K[i];
		}

		beta0 = aux1/(bn2*bn0-bn1*bn1) + beta0;
		beta1 = aux2/(bn2*bn0-bn1*bn1) + beta1;

		// Dispersion part
		for(int i=0; i<(n-1); ++i){
			tsG(i) = beta0G+beta1G*s_kk(i);
			etsG(i) = exp(tsG(i));
		}
		
		gi_old = gi;
		gi = my_fun_PoisM(x,x,y,gi_old,kappa,exp(psii1));
		gi0_kk = gi[0];
		gi0_kk.erase(kk);
		di = my_fun_d(y_kk,gi0_kk);
		q = tsG + di/etsG - 1;
	
		sn0G = 0;
		sn1G = 0;
		sn2G = 0;
		for(int i=0; i<(n-1); ++i){
			sn0G += kkxtG[i];
			sn1G += kkxtG[i]*s_kk[i];
			sn2G += kkxtG[i]*s_kk[i]*s_kk[i];
		}

		for(int i=0; i<(n-1); ++i){
			bjti0G(i) = kkxtG(i)*(sn2G-s_kk(i)*sn1G);
			bjti1G(i) = kkxtG(i)*(s_kk(i)*sn0G-sn1G);		
		}
	
		aux1 = 0;
		aux2 = 0;
		for(int i=0; i<(n-1); ++i){
			aux1 += bjti0G[i]*q[i];
			aux2 += bjti1G[i]*q[i];	
		}
		beta0G = aux1/(sn2G*sn0G-sn1G*sn1G);
		beta1G = aux2/(sn2G*sn0G-sn1G*sn1G);

		resi = {abs(beta0G/beta0G_old),abs(beta1G/beta1G_old),abs(beta0/beta0_old),abs(beta1/beta1_old)};
		res = max(resi);
		numit = numit + 1;

	}

	List L;
	if(numit == maxit  ){
		L = List::create(1);
	}else{
		L = List::create(beta0,beta0G);
	}
	

	return L;
	
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double my_fun_loglik_cv_nu(NumericVector x, NumericVector y, NumericVector startv, NumericVector startv2, const double kappa, const double nu){
	
	int n = y.size();
	NumericVector error(n);
	
	
	List fit;
	double fitM;
	double fitG;
	double di;
	int bien = 1;
	int m;
	for(int i=0; i<n; ++i){
		
		fit = my_fun_DoublePois_j_v2(x[i],x,y,startv,startv2,kappa,nu,i+1);
		m = fit.size();
		if(m==1){
			bien = 2;
			break;
		}else{
			
			fitM = fit[0];
			fitG = fit[1];
			if(isnan(fitM) || isnan(fitG)){
				bien = 2;
				break;	
			}else{
				di = my_fun_di(y(i),fitM);
				error[i] = fitG + di/exp(fitG);	
			}
		}
		
		
	
	}
	
	double res;
	if(bien==1){
		res = -sum(error)/2;
	}else{
		res = std::numeric_limits<double>::quiet_NaN();
	}
	
	
	return(res);

}
	







