function [H,gamma_o] = NL_RLDA(X,mu,n_tr)

% The function returns an estimate of the inverse covariance matrix that
% minimizes the misclassification probability.
%
% Input:
%       X:    is the raw data matrix of size n x p
%             n is the number of training samples.
%             p is the dimesionality.
%       mu:    matrix with class means of training data set
%              of size p x K) as column vectors
%       n_tr:  A vector of number of samples in each class. 

% Outputs:           
%          H        : estimate of the inverse covariance matrix   
%          gamma_o  : optimum gamma
%
% Reference: 
%  Maaz Mahadi "Regularized Linear Discriminant Analysis Using a
%  Nonlinear Covariance Matrix Estimator"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,p]    = size(X);    % data matrix dimensions

n_tilde   = n-2;

S         = X'*X/(n-2); % Sample Covariance Matrix

mu_hat    = mu';

m0        = mu_hat(:,1)   ;

m1        = mu_hat(:,2)   ;

m_plus    = m0  + m1;

m_minus   = m0  - m1;

n0        = n_tr(1);

n1        = n_tr(2);

c0        = n_tr(1)/n;       %prior probability of class 0
    
c1        = n_tr(2)/n;       %prior probability of class 1
    
tau       = log(c1/c0);    
    
gamma_vec = 1e5.^(linspace(-10,10,21)./10);   % range values of gamma
    
P_err_vec = nan(1,length(gamma_vec));         % intialize a vecotr estimating the error probability at each gamma 

if (p < n)   

    for kk  = 1:length(gamma_vec)
            H       =    1/gamma_vec(kk)^2*(inv(eye(p) + 1/gamma_vec(kk)*S))^2*S;
            Q     =    inv(gamma_vec(kk)*eye(p) + S);            
            G0g     =    (m0 - 1/2*m_plus)'*H*m_minus;
            G1g     =    (m1 - 1/2*m_plus)'*H*m_minus;
            Dg11    =    m_minus'*Q*S*Q*m_minus;
            Dg22    =    m_minus'*Q^2*S*Q^2*m_minus;    
            Dg12    =    m_minus'*Q^2*S*Q*m_minus;            
            [em_hat,em_hat_p1,theta1_hat] = CE(n,p, gamma_vec(kk),S);  % Consistent Estimator.           
            T2      =    em_hat_p1*(1+em_hat)*Dg11 + (1+em_hat)^2*Dg12;
            Dc    =    gamma_vec(kk)^2*(1+em_hat)^4*Dg22-2*gamma_vec(kk)*T2+(1+em_hat)^2*Dg11;
            if (Dc>=0) 
                sqrtDc = sqrt(Dc);
            else
                continue  % avoid the negative values of Dc
            end                 
            G0             =  (-G0g  + 1/n0*theta1_hat);
            G1             =  ( G1g  + 1/n1*theta1_hat);                      
            epslon_0       =  normcdf((G0+tau)/sqrtDc,0,1);  %e0
            epslon_1       =  normcdf((G1-tau)/sqrtDc,0,1);  %e1
            P_err_vec(kk)  =  c0*epslon_0 + c1*epslon_1;   % estimated total error probability 
    end
    [~,indx_opt]     =  min(P_err_vec);
    gamma_o          =  gamma_vec(indx_opt);
    H                =  1/gamma_o^2*(inv(eye(p) + 1/gamma_o*S))^2*S;
else
           [~,D,V] = svd(X,'econ');
    for kk = 1:length(gamma_vec)
            Q       = (D^2/n_tilde + gamma_vec(kk)*eye(n))^(-1);   
            H        = V*Q^2*D^2*V'/n_tilde;   
            G0g      = (m0 - 1/2*m_plus)'*H*m_minus;
            G1g      = (m1 - 1/2*m_plus)'*H*m_minus;
            H22      = V*Q^4*D^2*V'/n_tilde;
            H12      = V*Q^3*D^2*V'/n_tilde;
            Dg11     = m_minus'*H*m_minus;
            Dg22     = m_minus'*H22*m_minus;    
            Dg12     = m_minus'*H12*m_minus;            
            [em_hat,em_hat_p1,theta1_hat] = CE(n,p,gamma_vec(kk),D); % Consistent Estimator.              
            T2       = em_hat_p1*(1+em_hat)*Dg11 + (1+em_hat)^2*Dg12;
            Dc       = gamma_vec(kk)^2*(1+em_hat)^4*Dg22-2*gamma_vec(kk)*T2+(1+em_hat)^2*Dg11;
            if (Dc>=0) 
                sqrtDc = sqrt(Dc);
            else
                continue % avoid the negative values of Dc
            end              
            G0             =  (-G0g  + 1/n0*theta1_hat);
            G1             =  ( G1g  + 1/n1*theta1_hat);                      
            epslon_0       =  normcdf((G0+tau)/sqrtDc,0,1);  %e0
            epslon_1       =  normcdf((G1-tau)/sqrtDc,0,1);  %e1
            P_err_vec(kk)  =  c0*epslon_0 + c1*epslon_1;   %total error (Proposed consistent)
    end
    [~,indx_opt]       =  min(P_err_vec);
    gamma_o            =  gamma_vec(indx_opt);
     H                 =  V*(D^2/n_tilde + gamma_o*eye(n))^(-2)*D^2*V'/n_tilde;
end