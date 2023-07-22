function [em_hat,em_hat_p1b, theta_G_hat] = CECEP_LC(n0,n1,p,gama,D)
%clear all

n=n0+n1;
z = -gama;
I_n = eye(n);
  n=n-2;
    Hb = (D^2/n - z*I_n)^(-1);
    fzp1   = 1/n*trace ( Hb*D^2/n); 
    fzp2   = 1/n*trace ( Hb^2*D^2/n);
    em_hat     = real(fzp1/(1-fzp1));
    em_hat_p1b = real((1+em_hat)^2*fzp2 );
    xm_hat     =  1/(1+em_hat);
    xm_hat_p1  =  -em_hat_p1b/(1+em_hat)^2;
    theta_G_hat   = n*(em_hat_p1b*(xm_hat-z*xm_hat_p1) -fzp2 )/(em_hat_p1b*(1+em_hat)^(-2));







