function [em_hat,em_hat_p1, theta_G_hat] = CECEP(n,p,gamma,S)
% The function returns three estimated qunatities: 
% Output:
%        em_hat     :  Eq.(34) in the paper.
%        em_hat_p1  :  derivative of em_hat.
%        theta_G_hat: Eq. (33)
%
% Input: 
%        n          : number of training samples.
%        p          : dimensionality.
%        S          : Sample covariance matrix (or matrix of singular values)
% Reference: 
%  Maaz Mahadi "Regularized Linear Discriminant Analysis Using a
%  Nonlinear Covariance Matrix Estimator"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z        = -gamma;
n_tilde  = n-2   ;

if (p < n)
    fzp1          =  1/n_tilde*trace ( S*(S - z*eye(p))^(-1) );
    
    fzp2          =  1/n_tilde*trace ( S*(S - z*eye(p))^(-2) );
else
    Hb            = (S^2/n_tilde - z*eye(n))^(-1);
    
    fzp1          = 1/n_tilde*trace ( Hb*S^2/n_tilde); 
    
    fzp2          = 1/n_tilde*trace ( Hb^2*S^2/n_tilde);
end

em_hat        = real(fzp1/(1-fzp1));
    
em_hat_p1     = real((1+em_hat)^2*fzp2 );
    
xm_hat        =  1/(1+em_hat);
    
xm_hat_p1     =  -em_hat_p1/(1+em_hat)^2;
    
theta_G_hat   = n_tilde*(em_hat_p1*(xm_hat-z*xm_hat_p1) -fzp2 )/(em_hat_p1*(1+em_hat)^(-2));







