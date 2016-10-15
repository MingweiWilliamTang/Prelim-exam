#include "mcmc.head.h"

arma::mat TwoStateTransRate(double k, double pi1){
  double k12=k*pi1;
  double k21=k*(1-pi1);
  arma::mat Q;
  Q.zeros(2,2);
  Q(0,0) = -k12;
  Q(0,1) =  k12;
  Q(1,0) =  k21;
  Q(1,1) = -k21;
  return(Q);
}


arma::mat Hmatrix(double A0,double a, double b){
  arma::mat H;
  H.eye(2,2);
  H(0,0)=A0/a;
  H(1,1)=A0/b;
  return(H);
}


arma::mat Dmatrix(double tau,double a,double b){
  arma::mat D;
  D.eye(2,2);
  D(0,0)=a*exp(-a*tau);
  D(1,1)=b*exp(-b*tau);
  return(D);
}
//[[Rcpp::export()]]
arma::mat expm1(arma::mat A)
{
  arma::mat B;
  B.zeros(2,3);
  double t= (A(0,0)+A(1,1))/2;
  double d = sqrt((A(0,0)-A(1,1))*(A(0,0)-A(1,1))+4*A(0,1)*A(1,0));
  if (d>500)
  {
   B(0,0) =  (d*(1.0/2.0)+(A(0,0)-A(1,1))*(1.0/2.0))/d;
   B(0,1) = 2*A(0,1)*(1.0/2.0)/d;
   B(1,0) = 2*A(1,0)*(1.0/2)/d;
   B(1,1) = (d*(1.0/2)-(A(0,0)-A(1,1))*(1.0/2))/d;
   B(0,2) =t+d/2;
  }
  if( d<=500 && d!=0)
  {
  B(0,0) = (d*cosh(d/2)+(A(0,0)-A(1,1))*sinh(d/2))/d;
  B(0,1) = 2*A(0,1)*sinh(d/2)/d;
  B(1,0) = 2*A(1,0)*sinh(d/2)/d;
  B(1,1) = (d*cosh(d/2)-(A(0,0)-A(1,1))*sinh(d/2))/d;
  B(0,2) =t;
  }
  if(d==0)
  {
    B(0,0) = (1+(A(0,0)-A(1,1))/2);
    B(0,1) = A(0,1);
    B(1,0) = A(1,0);
    B(1,1) = (1-(A(0,0)-A(1,1))/2) ;
    B(0,2) =t;
  }
  return(B);
}


/*double ss(arma::vec a){
  return(a(2));
}*/



//the likelihood function for a single series a data
//[[Rcpp::export()]]

double likelihood_1(arma::vec dt,arma::vec tau,double a,double b,double pi1,double k, double A0)
{
  // scale the matrix take log
  double scale=0,s,t1; 
  arma::mat pis,tempmx;
  //tempmx stores the matrix temporaryily
  pis.zeros(1,2);
  pis(0,0) = pi1;
  pis(0,1) = 1-pi1;
  
  int n;
  n= dt.n_rows;
//  std::cout<<n<<endl;
  
  arma::mat Q = TwoStateTransRate(k,pi1);
  arma::mat H = Hmatrix(A0,a,b);
  arma::mat temp; // exp(Q-H)dt D_i+1 H
  temp.eye(2,2);
  arma::mat D,mQ,expmres;
  for(int i=0;i<=(n-1);i++)
  {
    D = Dmatrix(tau(i+1),a,b);
    expmres=expm1((Q-H)*dt(i));

        t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
     // Rcout<<mQ<<'\n'<<expmres<<endl;
    tempmx= temp*mQ*D*H; 
    s=fabs(tempmx(0,0))+fabs(tempmx(0,1))+fabs(tempmx(1,0))+fabs(tempmx(1,1));
    temp = tempmx/s;
    scale+=(log(s)+t1);
  }
  double L;
  arma::mat one;
  one.ones(2,1);
  arma::mat l = (pis*Dmatrix(tau(0),a,b)*temp*one);
  L=log(l(0,0))+scale;
  return(L);
}








//[[Rcpp::export()]]


arma::mat SimpleMCMC(arma::vec dt, arma::vec tau, arma::vec c,int iter, int burnin, int thin)
{ 
  //dt length vector
  //dt
  //cs is the vector associated with the parameter of the proposal distribution
  //total number of iteration = burnin + iter*thin
  //burnin  
  //cs =(c1,c2,c3,c4,c5)
  
  double a0,b0,pi10,k0,A00;
  double L0,L1,r;
  double temp1,temp2;
  arma::mat mcmcout(iter,6);
  
  int th;
  // initialize
  a0 = R::rgamma(1.0,2.0);
  b0 = R::rgamma(1.5625,1.0/1.5625);
  k0 = R::rexp(40000.0);
  pi10 = R::rbeta(0.89,0.89);
  A00 = R::rgamma(1.96,100000/5.6);
  double a=a0,b=b0,pi1=pi10,k=k0,A0=A00;
  // set a0>b)
a = (a0 > b0)?a0:b0;
b = (a0 <= b0)?a0:b0;
a0=a;
b0=b;

L0 = likelihood_1(dt,tau,a,b,pi1,k,A0);
int id;
for(int i=0; i<(iter+burnin); i++){
  th = ( (i<burnin) ? 1 : thin);
  for(int j=0;j<th;j++)
  {
    // sampling a b
    
    temp1 = R::rgamma(1/c(0),(a0*c(0)) ); 
    temp2 = R::rgamma(1/c(1),(b0*c(1)) );
    if(temp1>temp2) {
      a=temp1;
      b=temp2;
      id=1;
    }
    else {
      a=temp2;
      b=temp1;
      id=0;
    }  
    
    L1 = likelihood_1(dt,tau,a,b,pi10,k0,A00);
    if(id==1){
      r = L1+R::dgamma(a0,1.0/c(0),a*c(0),1)+R::dgamma(b0,1.0/c(1),b*c(1),1) +\
      R::dgamma(a,1.0,2.0,1)+R::dgamma(b,1.5625,1.0/1.5625,1)-\
      (L0 + R::dgamma(a,1.0/c(0),(a0*c(0)),1)+R::dgamma(b,1.0/c(1),(b0*c(1)),1)+\
       R::dgamma(a0,1.0,2.0,1)+R::dgamma(b0,1.5625,1.0/1.5625,1));
    }
    else
    {
      r = L1+R::dgamma(a0,1.0/c(0),b*c(0),1)+R::dgamma(b0,1.0/c(1),a*c(1),1) +\
      R::dgamma(a,1.0,2.0,1)+R::dgamma(b,1.5625,1.0/1.5625,1)-\
      (L0 + R::dgamma(b,1.0/c(0),(a0*c(0)),1)+R::dgamma(a,1.0/c(1),(b0*c(1)),1)+\
       R::dgamma(a0,1.0,2.0,1)+R::dgamma(b0,1.5625,1.0/1.5625,1));
    }
    r=exp(r);
    if(R::runif(0,1)<r) { a0=a;b0=b;  L0 = L1;}
    
    //sampling pi1
    pi1 = R::rbeta(c(2)*pi10,c(2)*(1.0-pi10));
    
    L1 = likelihood_1(dt,tau,a0,b0,pi1,k0,A00);
    
    r = L1 + R::dbeta(pi10,c(2)*pi1,c(2)*(1.0-pi1),1) + R::dbeta(pi1,0.89,0.89,1)-\
    (L0 + R::dbeta(pi1,c(2)*pi10,c(2)*(1.0-pi10),1) + R::dbeta(pi10,0.89,0.89,1));
    
    r=exp(r);
    if(R::runif(0,1)<r) { pi10 = pi1; L0 = L1;}
    
    
    //sampling k
    k = R::rgamma(1/c(3),(k0*c(3)));
    
    L1 = likelihood_1(dt,tau,a0,b0,pi10,k,A00);
    r = L1 + R::dgamma(k0,1/c(3),(k*c(3)),1) + R::dexp(k,40000.0,1) -\
    (L0 + R::dgamma(k,1/c(3),(k0*c(3)),1) + R::dexp(k0,40000.0,1) );
    r=exp(r);
    
    if(R::runif(0,1)<r) { k0 = k; L0 = L1;}
    
    // A0 
    A0 = R::rgamma(1/c(4),A00*c(4));
    L1 = likelihood_1(dt,tau,a0,b0,pi10,k0,A0);
    
    r = L1 + R::dgamma(A00,1/c(4),(A0*c(4)),1) + R::dgamma(A0,1.96,100000.0/5.6,1) -\
    (L0 + R::dgamma(A0,1/c(4),(A00*c(4)),1) + R::dgamma(A00,1.96,100000.0/5.6,1) );
    r=exp(r);
    if(R::runif(0,1)<r) { A00 = A0;  L0 = L1;}
  }
  //Rcout<<L0<<endl;
  //Rcout<<L1<<endl;
  //Rcout << a << '\t' <<b<<'\t'<< pi1 << '\t' << A0<<'\t'<<k<<endl;
  if(i>=burnin)
  {
    mcmcout(i-burnin,0) = a0;
    mcmcout(i-burnin,1) = b0;
    mcmcout(i-burnin,2) = pi10;
    mcmcout(i-burnin,3) = k0;
    mcmcout(i-burnin,4) = A00;
    mcmcout(i-burnin,5) = L0;
  }  
  
}
return(mcmcout);
}












//[[Rcpp::export()]]

double likelihood_2(arma::vec dt,arma::vec tau, arma::vec alpha, double a,double b,double pi1,double k, double A0){
  // scale the matrix take log
  double scale=0,s,t1; 
  arma::mat pis,tempmx;
  //tempmx stores the matrix temporaryily
  pis.zeros(1,2);
  pis(0,0) = pi1;
  pis(0,1) = 1-pi1;
  
  int n;
  n= dt.n_rows;
  //  std::cout<<n<<endl;
  
  arma::mat Q = TwoStateTransRate(k,pi1);
  arma::mat H = Hmatrix(A0,a,b);
  arma::mat Hi0(2,2),Hi1(2,2),mQ,expmres;
  arma::mat temp; // exp(Q-H)dt D_i+1 H
  temp.eye(2,2);
  arma::mat D;
  Hi0=alpha(0) * H;
  for(int i=0;i<=(n-1);i++)
  {
    Hi1 = alpha(i+1) * H;
    D = Dmatrix(tau(i+1),a,b);
    expmres=expm1((Q-Hi0)*dt(i));
    
    t1=expmres(0,2);
    mQ=expmres.submat(0,0,1,1); 
    tempmx= temp*mQ*D*Hi1; 
    s=fabs(tempmx(0,0))+fabs(tempmx(0,1))+fabs(tempmx(1,0))+fabs(tempmx(1,1));
    temp = tempmx/s;
    scale+=(log(s)+t1);
    Hi0 = Hi1;
  }
  double L;
  arma::mat one;
  one.ones(2,1);
  arma::mat l = (pis*Dmatrix(tau(0),a,b)*temp*one);
  L=log(l(0,0))+scale;
  return(L);
}

//[[Rcpp::export()]]

arma::vec thetaupdate1(arma::vec dt, arma::vec tau, arma::vec alpha, arma::vec c, double a0,double b0,double pi10, double k0, double A00,int iter)
{
  double a,b,k,pi1,A0,L1,temp1,temp2,r;
  double L0 = likelihood_2(dt,tau,alpha,a0,b0,pi10,k0,A00);
  arma::vec mcmcout(6);
  for(int i=0;i<iter;i++)
  {
    // sampling a b
    temp1 = R::rgamma(1/c(0),(a0*c(0)) ); 
    temp2 = R::rgamma(1/c(1),(b0*c(1)) );
    if(temp1>temp2) {
      a=temp1;
      b=temp2;
    }
    else {
      a=temp2;
      b=temp1; 
    }  
    L1 = likelihood_2(dt,tau,alpha,a,b,pi10,k0,A00);
    r = L1+R::dgamma(a0,1.0/c(0),(a*c(0)),1)+R::dgamma(b0,1.0/c(1),(b*c(1)),1) +\
    R::dgamma(a,1.0,2.0,1)+R::dgamma(b,1.5625,1.0/1.5625,1)-\
    (L0 + R::dgamma(a,1.0/c(0),(a0*c(0)),1)+R::dgamma(b,1.0/c(1),(b0*c(1)),1)+\
     R::dgamma(a0,1.0,2.0,1)+R::dgamma(b0,1.5625,1.0/1.5625,1));
    
    r=exp(r);
    if(R::runif(0,1)<r) { a0=a;b0=b;  L0 = L1;}
    
    //sampling pi1
    pi1 = R::rbeta(c(2)*pi10,c(2)*(1.0-pi10));
    
    L1 = likelihood_2(dt,tau,alpha,a0,b0,pi1,k0,A00);
    
    r = L1 + R::dbeta(pi10,c(2)*pi1,c(2)*(1.0-pi1),1) + R::dbeta(pi1,0.89,0.89,1)-\
    (L0 + R::dbeta(pi1,c(2)*pi10,c(2)*(1.0-pi10),1) + R::dbeta(pi10,0.89,0.89,1));
    
    r=exp(r);
    if(R::runif(0,1)<r) { pi10 = pi1; L0 = L1;}
    
    
    //sampling k
    k = R::rgamma(1/c(3),(k0*c(3)));
    
    L1 = likelihood_2(dt,tau,alpha,a0,b0,pi10,k,A00);
    r = L1 + R::dgamma(k0,1/c(3),(k*c(3)),1) + R::dexp(k,40000.0,1) -\
    (L0 + R::dgamma(k,1/c(3),(k0*c(3)),1) + R::dexp(k0,40000.0,1) );
    r=exp(r);
    
    if(R::runif(0,1)<r) { k0 = k; L0 = L1;}
    
    // A0 
    A0 = R::rgamma(1/c(4),A00*c(4));
    L1 = likelihood_2(dt,tau,alpha,a0,b0,pi10,k0,A0);
    
    r = L1 + R::dgamma(A00,1/c(4),(A0*c(4)),1) + R::dgamma(A0,1.96,100000.0/5.6,1) -\
    (L0 + R::dgamma(A0,1/c(4),(A00*c(4)),1) + R::dgamma(A00,1.96,100000.0/5.6,1) );
    r=exp(r);
    if(R::runif(0,1)<r) { A00 = A0;  L0 = L1;}
  }
  
  mcmcout(0) = a0;
  mcmcout(1) = b0;
  mcmcout(2) = pi10;
  mcmcout(3) = k0;
  mcmcout(4) = A00;
  mcmcout(5) = L0;
  return(mcmcout);
}  







//[[Rcpp::export()]]

double Alpha(double Bx,double By, double Bz,arma::vec w)
{
  return(exp(-((Bx*Bx)+(By*By))/(2*w(0)*w(0)) - (Bz*Bz)/(2*w(1)*w(1))) );
}

//[[Rcpp::export()]]
double Bratio(int i,arma::vec Bx,arma::vec dt,double Bx1,double sigma=1.0)
{
  int n = dt.n_rows;
  if(i==0) 
  {
    return(-1/(2*sigma*sigma*dt(i)) * ((Bx(i+1)-Bx1)*(Bx(i+1)-Bx1) -\
                                       (Bx(i+1)-Bx(i))*(Bx(i+1)-Bx(i))) );
  }
  else if (i == n)
  {
    return(-1/(2*sigma*sigma*dt(i-1)) * ((Bx(i-1)-Bx1)*(Bx(i-1)-Bx1) -\
                                         (Bx(i-1)-Bx(i))*(Bx(i-1)-Bx(i)) ) );   
  }
  else
  {
    return(-1/(2*sigma*sigma*dt(i)) * ( (Bx(i+1)-Bx1)*(Bx(i+1)-Bx1) -\
                                        (Bx(i+1)-Bx(i)) *(Bx(i+1)-Bx(i))) -\
           1/(2*sigma*sigma*dt(i-1)) * ((Bx(i-1)-Bx1)*(Bx(i-1)-Bx1) -\
                                        (Bx(i-1)-Bx(i)) * (Bx(i-1)-Bx(i)) ) ); 
  }
}


// [[Rcpp::export]]
List forbackwardMCMC(arma::vec dt, arma::vec tau, 
                     arma::vec c,arma::vec w,
                     double sigma, int iter, int burnin, int thin,
                     arma::vec BX,arma::vec BY,arma::vec BZ)
{
  int n = dt.n_rows;
  arma::mat paths(iter,n+1);
  double a0,b0,pi10,k0,A00;
  double a,b,pi1,k,A0;
  double Bx1,By1,Bz1,alpha1;
  double L0,L1,r;
  arma::mat expmres;
  double temp1,temp2,t1,t2;
  arma::vec Bx(n+1),By(n+1),Bz(n+1),alpha(n+1);
  arma::mat mcmcout(iter,6),D,mQ;
  arma::vec theta,scale2(n),scale(n);
  int th;
  
  //initialize 
  // theta  
  a0 = R::rgamma(1.0,2.0);
  b0 = R::rgamma(1.5625,1.0/1.5625);
  k0 = R::rexp(40000.0);
  pi10 = R::rbeta(0.89,0.89);
  A00 = R::rgamma(1.96,100000/5.6);
  a=a0,b=b0,pi1=pi10,k=k0,A0=A00;
  
  a = (a0 > b0)?a0:b0;
  b = (a0 <= b0)?a0:b0;
  a0=a;
  b0=b;
  
  // Bx,By,Bz
  Bx(0) = By(0) = Bz(0) = 0;
  
/*  {
    Bx=BX;
    By=BY;
    Bz=BZ;
  }*/
  
  alpha(0) = exp(-( (Bx(0)*Bx(0) +By(0)*By(0) )/ (2*w(0)*w(0))) - \
  ((Bz(0)*Bz(0))/(2*w(1)*w(1))));

  
  for(int i=1;i<n+1;i++)
  {
    Bx(i)=R::rnorm(Bx(i-1),sigma*sqrt(dt(i-1)));
    By(i)=R::rnorm(By(i-1),sigma*sqrt(dt(i-1)));
    Bz(i)=R::rnorm(Bz(i-1),sigma*sqrt(dt(i-1)));
    alpha(i)= exp(-( (Bx(i)*Bx(i) +By(i)*By(i) )/ (2*w(0)*w(0))) - \
                  ((Bz(i)*Bz(i))/(2*w(1)*w(1))));
  }
  // MCMC procedure
  for(int it=0; it<iter+burnin; it++)
  {
    if(it%500==0){Rcout<<it<<endl;}
    th = ( (it<burnin) ? 1 : thin);
    for(int j=0; j<th; j++)
    {
      // updata theta=c(a,b,pi,k,A0)
      
      theta=thetaupdate1(dt,tau,alpha,c,a0,b0,pi10,k0,A00,10);
      a = theta(0);
      b = theta(1);
      pi1 = theta(2);
      k = theta(3);
      A0 = theta(4);
      a0=a;b0=b;k0=k;
      pi10=pi1;A00=A0;
      
      
      //update Bx By Bz
      
      arma::mat pis;
      pis.zeros(1,2);
      pis(0,0) = pi1;
      pis(0,1) = 1-pi1;
      arma::mat Q = TwoStateTransRate(k,pi1);
      arma::mat H = Hmatrix(A0,a,b);
      
      // forward-backward
      // backward
      arma::cube Kcube(2,2,n+2);
      arma::vec scale(n+1),scale2(n+1);
      arma::mat Hi,Hi_til,temp,temp2; // temporarily store matrix Hi
      Kcube(0,0,n+1)=1;
      Kcube(1,0,n+1)=0;
      Kcube(0,1,n+1)=0;
      Kcube(1,1,n+1)=1;
      
      Kcube.slice(n) = alpha(n)*\
      Dmatrix(tau(n),a,b) * H;
      
      for(int i=n-1;i>=0;i--)
      {
        D = Dmatrix(tau(i),a,b);
        Hi = alpha(i) * H;
        expmres=expm1((Q-Hi) * dt(i));
        mQ = expmres.submat(0,0,1,1);
        
        temp = D * Hi * mQ * Kcube.slice(i+1);
        scale(i) = fabs(temp(0,0))+fabs(temp(1,1))+fabs(temp(0,1))+fabs(temp(1,0));
        Kcube.slice(i) = temp/scale(i);
      }
      
      //forward
      
      arma::mat R;
      arma::mat S;
      arma::mat vi=pis,one;
      arma::mat Lik;
      double l0,l1;
      one.ones(2,1);
      for(int i=0;i<=n;i++)
      {
        // proposed function for transition for Bx(i) to Bx1
        
        //Bx1 =R::rnorm(Bx(i),sigma*sqrt(d(i-1)));
        //By1 =R::rnorm(By(i),sigma*sqrt(d(i-1)));
        //Bz1 =R::rnorm(Bz(i),sigma*sqrt(d(i-1)));
        if(i==0) continue;
        else{
          if(it<2000){
          Bx1=R::rnorm(Bx(i-1),sigma*sqrt(dt(i-1)));
          By1=R::rnorm(By(i-1),sigma*sqrt(dt(i-1)));
          Bz1=R::rnorm(Bz(i-1),sigma*sqrt(dt(i-1)));
          }
          else{
            Bx1 =R::rnorm(Bx(i),c(5)*sigma*sqrt(dt(i-1)));
            By1 =R::rnorm(By(i),c(5)*sigma*sqrt(dt(i-1)));
            Bz1 =R::rnorm(Bz(i),c(5)*sigma*sqrt(dt(i-1)));
          }
          
          alpha1 = Alpha(Bx1,By1,Bz1,w);
          if(i<n)
          {
            D = Dmatrix(tau(i),a,b);
            Hi = alpha(i) * H;
            expmres=expm1((Q-Hi) * dt(i));
            mQ = expmres.submat(0,0,1,1);
            t1 = expmres(0,2);
            R = D * Hi* mQ;
            Hi_til =alpha1* H;
            expmres=expm1((Q-Hi_til)*dt(i)); 
            mQ = expmres.submat(0,0,1,1);
            t2 = expmres(0,2);
            S = D * Hi_til * mQ ;
            
          }
          else
          {
            D = Dmatrix(tau(i),a,b);
            Hi = alpha(i) * H;
            R = D * Hi;
            Hi_til = alpha1 * H;
            S = D * Hi_til;
            t1=t2=0;
          }
          // MH
          Lik = vi*R*Kcube.slice(i+1)*one;
          L0 = log(Lik(0,0))+t1;
          //    Rcout<<"L0"<<'\t'<<L0<<endl;
          Lik = vi*S*Kcube.slice(i+1)*one;
          L1 = log(Lik(0,0))+t2;
          //  Rcout<<"L1"<<'\t'<<L1<<endl;
          // Rcout<<L1-L0<<endl;
        //if(it>2000&& it<=4000){
        //r = L1-L0+(Bratio(i,Bx,dt,Bx1,sigma)+Bratio(i,By,dt,By1,sigma)+Bratio(i,Bz,dt,Bz1,sigma))/10;
         //}
         if(it>2000)  
        {r=L1-L0+(Bratio(i,Bx,dt,Bx1,sigma)+Bratio(i,By,dt,By1,sigma)+Bratio(i,Bz,dt,Bz1,sigma));}
        else{ r= L1-L0; }
        r = exp(r);
          //  Rcout<<Bx<<endl;
          //    Rcout<<"dif"<<alpha(i)-alpha1<<endl;
          //  Rcout<<i<<'\t'<<Bratio(i,Bx,dt,Bx1,sigma)<<endl;
          // Rcout<<"r"<<'\t'<<r<<endl;
          if(i==0){ r=-1;}
          if(R::runif(0,1) < r){
            Bx(i)=Bx1;
            By(i)=By1;
            Bz(i)=Bz1;
            temp2= vi*S;
            scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
            vi= temp2/scale2(i);
            alpha(i)=alpha1;
          }
          else{
            temp2= vi*R;
            scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
            vi= temp2/scale2(i); 
          }
        }
      }
    }
    if(it>=burnin){
      mcmcout(it-burnin,0) = a0;
      mcmcout(it-burnin,1) = b0;
      mcmcout(it-burnin,2) = pi10;
      mcmcout(it-burnin,3) = k0;
      mcmcout(it-burnin,4) = A00;
      mcmcout(it-burnin,5) = L0;
      paths.row(it-burnin)=alpha.t();
    }
  }
  List result;
  result["MCMC"]=mcmcout;
  result["path"]=paths;
  return(result);
}

arma::mat twostatesdf(double k, double pi1,double x)
{
  double k12=k*pi1;
  double k21=k*(1-pi1);
  arma::mat Q(2,2);
  Q(0,0) = -k12;
  Q(0,1) =  k12;
  Q(1,0) =  k21;
  Q(1,1) = -k21;
  Q=exp(-x) * Q;
  return(Q);
}

//[[Rcpp::export()]]

double likelihood_3(arma::vec dt,arma::vec tau, arma::vec alpha,arma::vec xt, 
                    double a,double b,double pi1,double k, double A0)
{
  
  // scale the matrix take log
  double scale=0,s,t1; 
  arma::mat pis,tempmx;
  //tempmx stores the matrix temporaryily
  pis.zeros(1,2);
  pis(0,0) = pi1;
  pis(0,1) = 1-pi1;
  int n;
  n= dt.n_rows;
  //  std::cout<<n<<endl;
  
  arma::mat Q(2,2);
  arma::mat H = Hmatrix(A0,a,b),expmres;
  arma::mat Hi0(2,2),Hi1(2,2),mQ;
  arma::mat temp; // exp(Q-H)dt D_i+1 H
  temp.eye(2,2);
  arma::mat D;
  Hi0=alpha(0) * H;
  for(int i=0;i<=(n-1);i++)
  {
    Hi1 = alpha(i+1) * H;
    D = Dmatrix(tau(i+1),a,b);
    Q = twostatesdf(k,pi1,xt(i));
    expmres=expm1((Q-Hi0)*dt(i));
    mQ = expmres.submat(0,0,1,1);
      t1 = expmres(0,2);
    tempmx= temp*mQ*D*Hi1; 
    s=fabs(tempmx(0,0))+fabs(tempmx(0,1))+fabs(tempmx(1,0))+fabs(tempmx(1,1));
    temp = tempmx/s;
    scale+=(log(s)+t1);
    Hi0 = Hi1;
  }
  double L;
  arma::mat one;
  one.ones(2,1);
  arma::mat l = (pis*Dmatrix(tau(0),a,b)*temp*one);
  L=log(l(0,0))+scale;
  return(L);
}



//[[Rcpp::export()]]

arma::vec thetaupdate2(arma::vec dt, arma::vec tau, arma::vec alpha, arma::vec xt, 
                       arma::vec c,
                       double a0,double b0,double pi10, double k0, double A00,
                       int iter)
{
  double a,b,k,pi1,A0,L1,temp1,temp2,r;
  double L0 = likelihood_3(dt,tau,alpha,xt,a0,b0,pi10,k0,A00);
  arma::vec mcmcout(6);
  for(int i=0;i<iter;i++)
  {
    // sampling a b
    temp1 = R::rgamma(1/c(0),(a0*c(0)) ); 
    temp2 = R::rgamma(1/c(1),(b0*c(1)) );
    if(temp1>temp2) {
      a=temp1;
      b=temp2;
    }
    else {
      a=temp2;
      b=temp1; 
    }  
    L1 = likelihood_3(dt,tau,alpha,xt,a,b,pi10,k0,A00);
    r = L1+R::dgamma(a0,1.0/c(0),(a*c(0)),1)+R::dgamma(b0,1.0/c(1),(b*c(1)),1) +\
    R::dgamma(a,1.0,2.0,1)+R::dgamma(b,1.5625,1.0/1.5625,1)-\
    (L0 + R::dgamma(a,1.0/c(0),(a0*c(0)),1)+R::dgamma(b,1.0/c(1),(b0*c(1)),1)+\
     R::dgamma(a0,1.0,2.0,1)+R::dgamma(b0,1.5625,1.0/1.5625,1));
    
    r=exp(r);
    if(R::runif(0,1)<r) { a0=a;b0=b;  L0 = L1;}
    
    //sampling pi1
    pi1 = R::rbeta(c(2)*pi10,c(2)*(1.0-pi10));
    
    L1 = likelihood_3(dt,tau,alpha,xt,a0,b0,pi1,k0,A00);
    
    r = L1 + R::dbeta(pi10,c(2)*pi1,c(2)*(1.0-pi1),1) + R::dbeta(pi1,0.89,0.89,1)-\
    (L0 + R::dbeta(pi1,c(2)*pi10,c(2)*(1.0-pi10),1) + R::dbeta(pi10,0.89,0.89,1));
    
    r=exp(r);
    if(R::runif(0,1)<r) { pi10 = pi1; L0 = L1;}
    
    
    //sampling k
    k = R::rgamma(1/c(3),(k0*c(3)));
    
    L1 = likelihood_3(dt,tau,alpha,xt,a0,b0,pi10,k,A00);
    r = L1 + R::dgamma(k0,1/c(3),(k*c(3)),1) + R::dexp(k,40000.0,1) -\
    (L0 + R::dgamma(k,1/c(3),(k0*c(3)),1) + R::dexp(k0,40000.0,1) );
    r=exp(r);
    
    if(R::runif(0,1)<r) { k0 = k; L0 = L1;}
    
    // A0 
    A0 = R::rgamma(1/c(4),A00*c(4));
    L1 = likelihood_3(dt,tau,alpha,xt,a0,b0,pi10,k0,A0);
    
    r = L1 + R::dgamma(A00,1/c(4),(A0*c(4)),1) + R::dgamma(A0,1.96,100000.0/5.6,1) -\
    (L0 + R::dgamma(A0,1/c(4),(A00*c(4)),1) + R::dgamma(A00,1.96,100000.0/5.6,1) );
    r=exp(r);
    if(R::runif(0,1)<r) { A00 = A0;  L0 = L1;}
  }
  
  mcmcout(0) = a0;
  mcmcout(1) = b0;
  mcmcout(2) = pi10;
  mcmcout(3) = k0;
  mcmcout(4) = A00;
  mcmcout(5) = L0;
  return(mcmcout);
}  

//[[Rcpp::export()]]
double logPath(arma::vec xt,arma::vec dt ,double lambda, double kxi)
{
  int n = xt.n_rows - 1;
  double s = -(n+1)/2 * log(kxi) - xt(0) * xt(0)/(2*kxi);
  for (int i=0;i<n;i++)
  {
    s += -log(1-exp(-2*lambda*dt(i)))/2 -\
    (xt(i+1)-xt(i)*exp(-lambda*dt(i)))*(xt(i+1)-xt(i)*exp(-lambda*dt(i)))/(2*kxi*(1-exp(-2*lambda*dt(i))));
  }
  return(s);
}




//[[Rcpp::export()]]
arma::vec UpdateOUPara(arma::vec xt, arma::vec dt, double lambda0, double kxi0, double c1, double c2,int iter = 10)
{
  arma::vec mcmcout(2);
  double lambda, kxi,r;
  double L0 = logPath(xt, dt ,lambda0 , kxi0), L1;
  for (int i=0; i<iter; i++)
  {
    lambda = R::rgamma(1/c1,c1* lambda0);
    L1 = logPath(xt, dt ,lambda, kxi0);
    r = L1 + R::dgamma(lambda0,1/c1,c1*lambda,1) + R::dgamma(lambda,40,5,1)-\
    L0- R::dgamma(lambda, 1/c1,c1* lambda0,1)-R::dgamma(lambda0,40,5,1);
    if (R::runif(0,1)< exp(r))
    {
      lambda0=lambda;
      L0=L1;
    }
    kxi = R::rgamma(1/c2,c2* kxi0);
    L1 = logPath(xt,dt, lambda0,kxi);
    r = L1 + R::dgamma(kxi0,1/c2,c2*kxi,1) + R::dgamma(kxi,1,2,1)-\
    L0- R::dgamma(kxi, 1/c2,c2* kxi0,1)-R::dgamma(kxi0,1,2,1);
    if ( R::runif(0,1)<exp(r))
    {
      kxi0=kxi;
      L0=L1;
    }
  }
  mcmcout(0)=lambda0;
  mcmcout(1)=kxi0;
  return(mcmcout);
}





//[[Rcpp::export()]]
double OUratio(int i, arma::vec xt,arma::vec dt, double xt1, double lambda, double kxi)
{
  double res;
  int n = dt.n_rows; 
  if(i==0)
  {
    res = -xt(0)*xt(0)/(2*kxi) -\
    (xt(i+1)-xt(i)*exp(-lambda*dt(i)))*(xt(i+1)-xt(i)*exp(-lambda*dt(i)))/(2*kxi*(1-exp(-2*lambda*dt(i))))+ \
    xt1*xt1/(2*kxi) +\
    (xt(i+1)-xt1*exp(-lambda*dt(i)))*(xt(i+1)-xt1*exp(-lambda*dt(i)))/(2*kxi*(1-exp(-2*lambda*dt(i))));
  }
  else if(i==n)
  {
    res = - (xt(i)-xt(i-1)*exp(-lambda*dt(i-1)))*(xt(i)-xt(i-1)*exp(-lambda*dt(i-1)))/(2*kxi*(1-exp(-2*lambda*dt(i-1)))) +\
    (xt1-xt(i-1)*exp(-lambda*dt(i-1)))*(xt1-xt(i-1)*exp(-lambda*dt(i-1)))/(2*kxi*(1-exp(-2*lambda*dt(i-1)))); \
  }  									
  else
  {
    res = - (xt(i)-xt(i-1)*exp(-lambda*dt(i-1)))*(xt(i)-xt(i-1)*exp(-lambda*dt(i-1)))/(2*kxi*(1-exp(-2*lambda*dt(i-1)))) - \
    (xt(i+1)-xt(i)*exp(-lambda*dt(i)))*(xt(i+1)-xt(i)*exp(-lambda*dt(i)))/(2*kxi*(1-exp(-2*lambda*dt(i))))+ \
    (xt1-xt(i-1)*exp(-lambda*dt(i-1)))*(xt1-xt(i-1)*exp(-lambda*dt(i-1)))/(2*kxi*(1-exp(-2*lambda*dt(i-1))))+ \
    (xt(i+1)-xt1*exp(-lambda*dt(i)))*(xt(i+1)-xt1*exp(-lambda*dt(i)))/(2*kxi*(1-exp(-2*lambda*dt(i))));
  }								     
  return(-res);
}






 //[[Rcpp::export()]]

List ContdiffMCMC(arma::vec dt, arma::vec tau, 
arma::vec c,arma::vec w,
double sigma, int iter, int burnin, int thin,
arma::vec BX,arma::vec BY,arma::vec BZ)
{
    int n = dt.n_rows;
    arma::mat paths(iter,n+1),xts(iter,n+1);
    double a0,b0,pi10,k0,A00,lambda0,kxi0;
    double a,b,pi1,k,A0,lambda,kxi;
    double Bx1,By1,Bz1,alpha1,xt1,t1,t2;
    double L0,L1,r;
    double temp1,temp2;
    arma::vec Bx(n+1),By(n+1),Bz(n+1),alpha(n+1),xt(n+1);
    arma::mat mcmcout(iter,8),D,expmres,mQ,Q_til;
    arma::vec theta,OUpara(2),scale2(n),scale(n),scale3(n),scale4(n);
  int th;

//initialize 
   // theta  
  a0 = R::rgamma(1.0,2.0);
  b0 = R::rgamma(1.5625,1.0/1.5625);
  k0 = R::rexp(40000.0);
  pi10 = R::rbeta(0.89,0.89);
  A00 = R::rgamma(1.96,100000/5.6);
   a=a0,b=b0,pi1=pi10,k=k0,A0=A00;
  
   a = (a0 > b0)?a0:b0;
   b = (a0 <= b0)?a0:b0;
   a0=a;
   b0=b;
  
  // Bx,By,Bz

Bx(0) = By(0) = Bz(0) = 0;
  alpha(0) = exp(-( (Bx(0)*Bx(0) +By(0)*By(0) )/ (2*w(0)*w(0))) - \
  ((Bz(0)*Bz(0))/(2*w(1)*w(1))));
  for(int i=1;i<n+1;i++)
  {
    Bx(i)=R::rnorm(Bx(i-1),sigma*sqrt(dt(i-1)));
    By(i)=R::rnorm(By(i-1),sigma*sqrt(dt(i-1)));
    Bz(i)=R::rnorm(Bz(i-1),sigma*sqrt(dt(i-1)));
    alpha(i)= exp(-( (Bx(i)*Bx(i) +By(i)*By(i) )/ (2*w(0)*w(0))) - \
  ((Bz(i)*Bz(i))/(2*w(1)*w(1))));
  }
  
/*  
  Bx=BX;
  By=BY;
  Bz=BZ;
 */ 
  alpha(0) = exp(-( (Bx(0)*Bx(0) +By(0)*By(0) )/ (2*w(0)*w(0))) - \
  ((Bz(0)*Bz(0))/(2*w(1)*w(1))));
  for(int i=1;i<n+1;i++)
  {
  alpha(i)= exp(-( (Bx(i)*Bx(i) +By(i)*By(i) )/ (2*w(0)*w(0))) - \
  ((Bz(i)*Bz(i))/(2*w(1)*w(1))));
  }
  
  // lambda xi
  lambda0 = R::rgamma(100 , 5.0);
  kxi0 = R::rgamma(1.0 , 2.0);
//  kxi=kxi0=0.1;
  lambda=lambda0;

  // xt
 // xt(0) = R::rnorm(0 , sqrt(kxi0));
  xt(0) = 0;
  for (int i = 1; i< n+1 ; i++)
    {
      xt(i) = R::rnorm(xt(i-1) * exp(-lambda * dt(i-1)),sqrt(kxi *\
      (1 - exp(-2*lambda * dt(i-1) ))));
    }
 
Rcout<<"mcmcstart"<<endl;

// MCMC procedure

  for(int it=0; it<iter+burnin; it++)
{
   th = ( (it<burnin) ? 1 : thin);
   if(it%500==0) {Rcout<<"finish "<<it<<" iterations"<<endl;}
   for(int j=0; j<th; j++)
   {
     
     

  // updata theta=c(a,b,pi,k,A0)
  
     theta=thetaupdate2(dt,tau,alpha,xt,c,a0,b0,pi10,k0,A00,20);
  a = theta(0);
  b = theta(1);
  pi1 = theta(2);
  k = theta(3);
  A0 = theta(4);
  a0=a;b0=b;k0=k;
  pi10=pi1; A00=A0;
  
  
  //updata kxi, lambda
   
  OUpara = UpdateOUPara(xt,dt,lambda0,kxi0,c(5),c(6));
  lambda0 = OUpara(0);
  kxi0 = OUpara(1);
  lambda =lambda0;
  kxi = kxi0;


  //update Bx By Bz



  arma::mat pis;
    pis.zeros(1,2);
  pis(0,0) = pi1;
  pis(0,1) = 1-pi1;
  arma::mat Q;
  arma::mat H = Hmatrix(A0,a,b);

  // forward-backward
  // backward
  arma::cube Kcube(2,2,n+2);
  arma::vec scale(n+1),scale2(n+1);
  arma::mat Hi,Hi_til,temp,temp2; // temporarily store matrix Hi
  Kcube(0,0,n+1)=1;
  Kcube(1,0,n+1)=0;
  Kcube(0,1,n+1)=0;
  Kcube(1,1,n+1)=1;
  
  Kcube.slice(n) = alpha(n)*\
  Dmatrix(tau(n),a,b) * H;
  
  for(int i=n-1;i>=0;i--)
  {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q = twostatesdf(k, pi1, xt(i));
    
    expmres=expm1((Q-Hi) * dt(i)) ;
        t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    
    temp = D * Hi* mQ * Kcube.slice(i+1);
    scale(i) = fabs(temp(0,0))+fabs(temp(1,1))+fabs(temp(0,1))+fabs(temp(1,0));
    Kcube.slice(i) = temp/scale(i);
  }
  
  //forward
  
  arma::mat R;
  arma::mat S;
  arma::mat vi = pis,one;
  arma::mat Lik;
  double l0,l1;
  one.ones(2,1);
  for(int i=0;i<=n;i++)
  {
    // proposed function for transition for Bx(i) to Bx1
    

    if(i==0) continue;
    else{
    Bx1=R::rnorm(Bx(i-1),sigma*sqrt(dt(i-1)));
    By1=R::rnorm(By(i-1),sigma*sqrt(dt(i-1)));
    Bz1=R::rnorm(Bz(i-1),sigma*sqrt(dt(i-1)));
/*
    Bx1=R::rnorm(Bx(i),sigma*sqrt(dt(i-1)));
    By1=R::rnorm(By(i),sigma*sqrt(dt(i-1)));
    Bz1=R::rnorm(Bz(i),sigma*sqrt(dt(i-1)));
*/  
    alpha1 = Alpha(Bx1,By1,Bz1,w);
   
    if(i<n)
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q = twostatesdf(k,pi1,xt(i));
    
    expmres=expm1((Q-Hi) * dt(i)) ;
    t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    
    R = D * Hi * mQ;
    Hi_til =alpha1  * H;
    
    expmres=expm1((Q-Hi_til) * dt(i)) ;
    t2=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    S = D * Hi_til * mQ;
    }
    else
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    R = D * Hi;
    Hi_til = alpha1 * H;
    S = D * Hi_til;
    }
// MH
    Lik = vi*R*Kcube.slice(i+1)*one;
    L0 = log(Lik(0,0))+t1;
//    Rcout<<"L0"<<'\t'<<L0<<endl;
    Lik = vi*S*Kcube.slice(i+1)*one;
    L1 = log(Lik(0,0))+t2;
  //  Rcout<<"L1"<<'\t'<<L1<<endl;
   // Rcout<<L1-L0<<endl;
    r = L1-L0;
    //+Bratio(i,Bx,dt,Bx1,sigma)+Bratio(i,By,dt,By1,sigma)+Bratio(i,Bz,dt,Bz1,sigma);
    r = exp(r);

    if(i==0){ r=-1;}
    if(R::runif(0,1) < r){
      Bx(i)=Bx1;
      By(i)=By1;
      Bz(i)=Bz1;
      temp2= vi*S;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i);
      alpha(i)=alpha1;
    }
    else{
      temp2= vi*R;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i); 
    }
    }
  }

  // updata xt

  //backward forward

 // backward
  Kcube(0,0,n+1)=1;
  Kcube(1,0,n+1)=0;
  Kcube(0,1,n+1)=0;
  Kcube(1,1,n+1)=1;
  vi = pis;
  Kcube.slice(n) = alpha(n)*\
  Dmatrix(tau(n),a,b) * H;
  
  for(int i=n-1;i>=0;i--)
  {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q = twostatesdf(k, pi1, xt(i));
    expmres=expm1((Q-Hi) * dt(i)) ;
    t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    temp = D * Hi * mQ * Kcube.slice(i+1);
    scale(i) = fabs(temp(0,0))+fabs(temp(1,1))+fabs(temp(0,1))+fabs(temp(1,0));
    Kcube.slice(i) = temp/scale(i);
  }
  

 //forward
  
 
  for(int i=0;i<=n;i++)
  {
    // proposed function for transition from xt(i) to xt1
    

    //xt1 = R::rgamma(1/c(7),c(7)*xt(i));
    xt1 = R::rnorm(xt(i),c(7));
    if(i<n)
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q = twostatesdf(k,pi1,xt(i));
    Q_til = twostatesdf(k,pi1,xt1);
    expmres=expm1((Q-Hi) * dt(i)) ;
    t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    R = D * Hi * mQ;
  //  Hi_til =alpha1  * H;
    expmres=expm1((Q_til-Hi) * dt(i)) ;
    t2=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    S = D * Hi * mQ;
    }
    else
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    R = D * Hi;
//    Hi_til = alpha1 * H;
    S = D * Hi;
    }
// MH
    Lik = vi*R*Kcube.slice(i+1)*one;
    L0 = log(Lik(0,0))+t1;
//    Rcout<<"L0"<<'\t'<<L0<<endl;
    Lik = vi*S*Kcube.slice(i+1)*one;
    L1 = log(Lik(0,0))+t2;
  //  Rcout<<"L1"<<'\t'<<L1<<endl;
    r = L1-L0+OUratio(i,xt,dt,xt1,lambda,kxi);\
    //R::dgamma(xt(i),1/c(7),c(7)*xt1,1)-R::dgamma(xt1,1/c(7),c(7)*xt(i),1);
    r = exp(r);

    if(R::runif(0,1) < r){
      xt(i)=xt1;
      temp2= vi*S;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i);
    }
    else{
      temp2= vi*R;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i); 
    }
   }
  }
   // store mcmc results

  if (it >= burnin) {
  mcmcout(it-burnin, 0) = a0;
  mcmcout(it-burnin, 1) = b0;
  mcmcout(it-burnin, 2) = pi10;
  mcmcout(it-burnin, 3) = k0;
  mcmcout(it-burnin, 4) = A00;
  mcmcout(it-burnin, 5) = lambda0;
  mcmcout(it-burnin, 6) = kxi0;
  mcmcout(it-burnin, 7) = likelihood_3(dt,tau,alpha,xt,a0,b0,pi10,k0,A00);
  paths.row(it-burnin) = alpha.t();
  xts.row(it-burnin) = xt.t();
}
}
List result;
result["MCMC"] = mcmcout;
result["path"] = paths;
result["energy"] = xts;
return(result);

}


  //[[Rcpp::export()]]
arma::mat Balpha(arma::vec dt,double sigma=1000)
{
  int n = dt.n_rows;
  arma::vec Bx(n+1),By(n+1),Bz(n+1),alpha(n+1),w(2);
  w(0) = 310; 
  w(1) = 1760;
  arma::mat res(4,n+1);
  Bx(0) = By(0) = Bz(0) = 0;
  alpha(0) = exp(-( (Bx(0)*Bx(0) +By(0)*By(0) )/ (2*w(0)*w(0))) - \
  ((Bz(0)*Bz(0))/(2*w(1)*w(1))));
  for(int i=1;i<n+1;i++)
  {
    Bx(i)=R::rnorm(Bx(i-1),sigma*sqrt(dt(i-1)));
    By(i)=R::rnorm(By(i-1),sigma*sqrt(dt(i-1)));
    Bz(i)=R::rnorm(Bz(i-1),sigma*sqrt(dt(i-1)));
    alpha(i)= exp(-( (Bx(i)*Bx(i) +By(i)*By(i) )/ (2*w(0)*w(0))) - \
  ((Bz(i)*Bz(i))/(2*w(1)*w(1))));
  }

 res.row(0)= Bx.t();
 res.row(1) =By.t();
 res.row(2) = Bz.t();
 res.row(3) = alpha.t();
  return(res);
}




  //[[Rcpp::export()]]


arma::mat updateBrowian(arma::vec dt, arma::vec tau,arma::vec alpha0,
arma::vec Bx,arma::vec By, arma::vec Bz, arma::vec xt,
double a, double b, double pi1,double k,double A0,double c1,int iter=1,int maxlik=0){
int n = dt.n_rows;
  arma::vec alpha=alpha0,w(2);
  w(0)=310;
  w(1)=1760;
 arma::mat pis(1,2);
    arma::mat alphas(iter,n+1);
    pis(0,0) = pi1;
  pis(0,1) = 1-pi1;
  arma::mat Q,D;
  arma::mat H = Hmatrix(A0,a,b);
    double sigma=1000,t1,t2;
 arma::mat result(4,n+1),expmres,mQ;
  arma::cube Kcube(2,2,n+2);
  arma::vec scale(n+1),scale2(n+1);
  arma::mat Hi,Hi_til,temp,temp2; // temporarily store matrix Hi
    double Bx1,By1,Bz1,r,L0,L1,alpha1;
for(int k=0; k<iter; k++)
{
  // forward-backward
  // backward
  Kcube(0,0,n+1)=1;
  Kcube(1,0,n+1)=0;
  Kcube(0,1,n+1)=0;
  Kcube(1,1,n+1)=1;

  Kcube.slice(n) = alpha(n)*\
  Dmatrix(tau(n),a,b) * H;
  
  for(int i=n-1;i>=0;i--)
  {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q = twostatesdf(k, pi1, xt(i));
    expmres=expm1((Q-Hi)*dt(i));
        t2=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    temp = D * Hi * mQ * Kcube.slice(i+1);
    scale(i) = fabs(temp(0,0))+fabs(temp(1,1))+fabs(temp(0,1))+fabs(temp(1,0));
    Kcube.slice(i) = temp/scale(i);
  }
  
  //forward

  arma::mat R;
  arma::mat S;
  arma::mat vi = pis,one;
  arma::mat Lik;
  double l0,l1;
  one.ones(2,1);
  for(int i=0;i<=n;i++)
  {
    // proposed function for transition for Bx(i) to Bx1
    
    if(i==0) {
      Bx1=0;//R::rnorm(Bx(0),c1*sigma*sqrt(dt(i))/10);
      By1=0;//R::rnorm(By(0),c1*sigma*sqrt(dt(i))/10);
      Bz1=0;//R::rnorm(Bz(0),c1*sigma*sqrt(dt(i))/10);
    }
  else{
  if(maxlik==0){
    Bx1=R::rnorm(Bx(i),c1*sigma*sqrt(dt(i-1)));
    By1=R::rnorm(By(i),c1*sigma*sqrt(dt(i-1)));
    Bz1=R::rnorm(Bz(i),c1*sigma*sqrt(dt(i-1)));
  }
  else{
      Bx1=R::rnorm(Bx(i-1),sigma*sqrt(dt(i-1)));
      By1=R::rnorm(By(i-1),sigma*sqrt(dt(i-1)));
      Bz1=R::rnorm(Bz(i-1),sigma*sqrt(dt(i-1)));
  }
  }
    alpha1 = Alpha(Bx1,By1,Bz1,w);
   
    if(i<n)
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q = twostatesdf(k,pi1,xt(i));
    expmres=expm1((Q-Hi)*dt(i));
        t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
     
    R = D * Hi * mQ;
    Hi_til =alpha1  * H;
    expmres=expm1((Q-Hi_til)*dt(i));
        t2=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    S = D * Hi_til * mQ;
    }
    else
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    R = D * Hi;
    Hi_til = alpha1 * H;
    S = D * Hi_til;
    t1=t2;
    }
// MH
    Lik = vi*R*Kcube.slice(i+1)*one;
    L0 = log(Lik(0,0))+t1;
    Lik = vi*S*Kcube.slice(i+1)*one;
    L1 = log(Lik(0,0))+t2;
  //  Rcout<<"L1"<<'\t'<<L1<<endl;
   // Rcout<<L1-L0<<endl;
    r = L1-L0;
//    +(Bratio(i,Bx,dt,Bx1,sigma)+Bratio(i,By,dt,By1,sigma)+Bratio(i,Bz,dt,Bz1,sigma));
    
    r = exp(r);

    if(R::runif(0,1) < r){
      Bx(i)=Bx1;
      By(i)=By1;
      Bz(i)=Bz1;
      temp2= vi*S;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i);
      alpha(i)=alpha1;
    }
    else{
      temp2= vi*R;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i); 
    }
  }

// alphas.row(k) = alpha.t();
}
  result.row(0) = Bx.t();
  result.row(1) = By.t();
  result.row(2) = Bz.t(); 
  result.row(3) = alpha.t();
return(result);
}





  //[[Rcpp::export()]]


arma::vec updatext(arma::vec dt, arma::vec tau,arma::vec alpha,
arma::vec xt0,double lambda,double kxi,
double a, double b, double pi1,double k,double A0,double c1,int iter=1)
{  // updata xt

int n = dt.n_rows;
  arma::vec xt=xt0,w(2);
  w(0)=310;
  w(1)=1760;
 arma::mat pis(1,2);
    arma::mat xts(iter,n+1);
    pis(0,0) = pi1;
  pis(0,1) = 1-pi1;
  arma::mat Q,D,Q_til;
  arma::mat H = Hmatrix(A0,a,b);
    double sigma=1000,t1,t2;
  arma::cube Kcube(2,2,n+2);
  arma::vec scale(n+1),scale2(n+1);
  arma::mat Hi,temp,temp2,vi,expmres,mQ; // temporarily store matrix Hi
    double r,L0,L1,xt1;
for(int k=0; k<iter; k++)
{ 

  //backward forward
 // backward
  Kcube(0,0,n+1)=1;
  Kcube(1,0,n+1)=0;
  Kcube(0,1,n+1)=0;
  Kcube(1,1,n+1)=1;
  vi = pis;
  Kcube.slice(n) = alpha(n)*\
  Dmatrix(tau(n),a,b) * H;
  
  for(int i=n-1;i>=0;i--)
  {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q = twostatesdf(k, pi1, xt(i));
    expmres=expm1((Q-Hi)*dt(i));
        t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    temp = D * Hi * mQ * Kcube.slice(i+1);
    scale(i) = fabs(temp(0,0))+fabs(temp(1,1))+fabs(temp(0,1))+fabs(temp(1,0));
    Kcube.slice(i) = temp/scale(i);
  }
  

 //forward
  arma::mat R;
  arma::mat S;
  arma::mat vi = pis,one;
  arma::mat Lik;
  double l0,l1;
  one.ones(2,1);
 
  for(int i=0;i<=n;i++)
  {
    // proposed function for transition from xt(i) to xt1
    

    xt1 = R::rgamma(1/c1,c1*xt(i));
    //if(i==0) {xt1 = R::rnorm(xt(i),c1);}
    //else {xt1=R::rnorm(xt(i-1) * exp(-lambda * dt(i-1)),sqrt(kxi *\
    //  (1 - exp(-2*lambda * dt(i-1) ))));}
    
    if(i<n)
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    Q_til = twostatesdf(k,pi1,xt1);
    Q = twostatesdf(k,pi1,xt(i));
    
    expmres=expm1((Q-Hi)*dt(i));
        t1=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    R = D * Hi * mQ;
    
    expmres=expm1((Q_til-Hi)*dt(i));
        t2=expmres(0,2);
     mQ=expmres.submat(0,0,1,1); 
    S = D * Hi * mQ;
    }
    else
    {
    D = Dmatrix(tau(i),a,b);
    Hi = alpha(i) * H;
    R = D * Hi;
    S = R;
    t1=t2=0;
    }
// MH
    Lik = vi*R*Kcube.slice(i+1)*one;
    L0 = log(Lik(0,0))+t1;

    Lik = vi*S*Kcube.slice(i+1)*one;
    L1 = log(Lik(0,0))+t2;
  //  Rcout<<"L1"<<'\t'<<L1<<endl;
    r = L1-L0+OUratio(i,xt,dt,xt1,lambda,kxi);
    //R::dgamma(xt(i),1/c(7),c(7)*xt1,1)-R::dgamma(xt1,1/c(7),c(7)*xt(i),1);
    r = exp(r);

    if(R::runif(0,1) < r){
      xt(i)=xt1;
      temp2= vi*S;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i);
    }
    else{
      temp2= vi*R;
      scale2(i)=fabs(temp2(0,0))+fabs(temp2(0,1));
      vi= temp2/scale2(i); 
    }
   }

//   xts.row(k) = xt.t();
  }
  return(xt);
  }
  
  
  //[[Rcpp::export()]]
arma::vec OUpath(arma::vec dt, double lambda,double kxi)
{
  int n= dt.n_rows;
  arma::vec xt(n+1);
  //xt(0) = R::rnorm(0 , sqrt(kxi));
  xt(0) = 0;
  for (int i = 1; i< n+1 ; i++)
    {
      xt(i) = R::rnorm(xt(i-1) * exp(-lambda * dt(i-1)),sqrt(kxi *\
      (1 - exp(-2*lambda * dt(i-1) ))));
    }
 return(xt);
}


//[[Rcpp::export()]]
double log_Br(arma::vec dt,arma::vec Bx,arma::vec By,
arma::vec Bz,double sigma=1000)
{
  int n= Bx.n_rows;
  double L=0;
  for(int i=0;i<n-1;i++)
  {
     L+=-1/(2*sigma*sigma*dt(i))*(Bx(i+1)-Bx(i))*(Bx(i+1)-Bx(i));
     L+=-1/(2*sigma*sigma*dt(i))*(By(i+1)-By(i))*(By(i+1)-By(i));
     L+=-1/(2*sigma*sigma*dt(i))*(Bz(i+1)-Bz(i))*(Bz(i+1)-Bz(i));
  }
  return(L);
}


