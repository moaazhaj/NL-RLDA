function [em_hat,em_hat_p1b, theta_G_hat] = CE_D2(n0,n1,p,gama,Sigma)


n=n0+n1;
M = p;
N = n;
R = Sigma ;




z = -1/gama;


  I_N = eye(N);
  I_n = eye(n);
  I_M = eye(M);
  I_p = eye(p);
  
  n=n-2;
  N=N-2;
 % C = (Y*Y'/N);
C = R;




fzp1   = 1/n*trace ( C*(C - z*eye(p))^(-1) ); 
fzp2   = 1/n*trace ( C*(C - z*eye(p))^(-2) );



f_gamaC2 = 1/n*trace ( C*(gama*C + eye(p))^(-2) ); 

f_gamaC_prime = - 1/n*trace ( C^2*(gama*C + eye(p))^(-2) );
f_zC_prime = 1/n*trace ( C*(C - z*eye(p))^(-2) ); 



em_hat     = real(fzp1/(1-fzp1));
em_hat_p1a = real(f_zC_prime/ (1-fzp1)^2);
em_hat_p1b = real((1+em_hat)^2*fzp2 );

xm_hat     =  1/(1+em_hat);
xm_hat_p1  =  -em_hat_p1b/(1+em_hat)^2;




epsi_hat_2_In = fzp2/(xm_hat-z*xm_hat_p1);
phi_hat       = real(em_hat_p1b-epsi_hat_2_In)/(em_hat_p1b*(1+em_hat)^(-2));
DE2_hat       = (xm_hat-z*xm_hat_p1)*phi_hat;
theta_G_hat   = n*(em_hat_p1b*(xm_hat-z*xm_hat_p1) -fzp2 )/(em_hat_p1b*(1+em_hat)^(-2));


















