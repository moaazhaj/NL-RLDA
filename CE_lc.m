function [em_hat,em_hat_p1b, theta_G_hat] = CE_lc(n0,n1,p,gama,D)


n=n0+n1;
M = p;
N = n;


z = -1/gama;


  I_N = eye(N);
  I_n = eye(n);
  I_M = eye(M);
  I_p = eye(p);
  
  n=n-2;
  N=N-2;





Hb = (D^2/n - z*I_n)^(-1);

fzp1   = 1/n*trace ( Hb*D^2/n); 
fzp2   = 1/n*trace ( Hb^2*D^2/n);




em_hat     = real(fzp1/(1-fzp1));
em_hat_p1b = real((1+em_hat)^2*fzp2 );


xm_hat     =  1/(1+em_hat);
xm_hat_p1  =  -em_hat_p1b/(1+em_hat)^2;



epsi_hat_2_In = fzp2/(xm_hat-z*xm_hat_p1);
phi_hat       = real(em_hat_p1b-epsi_hat_2_In)/(em_hat_p1b*(1+em_hat)^(-2));
DE2_hat       = (xm_hat-z*xm_hat_p1)*phi_hat;
theta_G_hat   = n*(em_hat_p1b*(xm_hat-z*xm_hat_p1) -fzp2 )/(em_hat_p1b*(1+em_hat)^(-2));







