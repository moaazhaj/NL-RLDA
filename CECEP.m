function [em_hat,em_hat_p1b, theta_G_hat] = CECEP(n0,n1,p,gama,Sigma)
%clear all

n=n0+n1;
M = p;
N = n;
R = Sigma ;
z = - gama;
n=n-2;
N=N-2;
C = R;
fzp1   = 1/n*trace ( C*(C - z*eye(p))^(-1) ); 
fzp2   = 1/n*trace ( C*(C - z*eye(p))^(-2) );
em_hat     = real(fzp1/(1-fzp1));
em_hat_p1b = real((1+em_hat)^2*fzp2 );
xm_hat        =  1/(1+em_hat);
xm_hat_p1      =  -em_hat_p1b/(1+em_hat)^2;
theta_G_hat   = n*(em_hat_p1b*(xm_hat-z*xm_hat_p1) -fzp2 )/(em_hat_p1b*(1+em_hat)^(-2));







