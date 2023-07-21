function [HNL,gamma_o] = nlrlda_lc(X,mu,n_tr)
S_hat = X'*X/(sum(n_tr)-2);
[U,D,V] = svd(X,'econ');
[nt p]= size(X);
n = nt-2;
    mu_hat   = mu';
    cnt=0;
    m_plus   = mu_hat(:,1) + mu_hat(:,2);
    m_minus  = mu_hat(:,1) - mu_hat(:,2);  
    c0 = n_tr(1)/ (n_tr(1)+ n_tr(2));
    c1 = n_tr(2)/ (n_tr(1)+ n_tr(2));
    cc=log(c1/c0);    
     p = length(S_hat);
   % gg          = -10:0.25:-5;
    gg=linspace(-10,10,21);
    gamma = 1e5.^(gg./10);
    epsl_Pc = zeros(1,length(gamma));
           for kk = 1:length(gamma)
            Hb = (D^2/n + gamma(kk)*eye(nt))^(-1);
            HNL=V*Hb^2*D^2*V'/n;
            G0g=(mu_hat(:,1) - 1/2*m_plus)'*HNL*m_minus;
            G1g=(mu_hat(:,2) - 1/2*m_plus)'*HNL*m_minus;
            Dg = m_minus'*HNL*S_hat*HNL*m_minus;
            H22= V*Hb^4*D^2*V'/n;
            H12= V*Hb^3*D^2*V'/n;
            Dg11 = m_minus'*HNL*m_minus;
            Dg22 = m_minus'*H22*m_minus;    
            Dg12 = m_minus'*H12*m_minus;            
            [em_hat,em_hat_p1,theta1_hat] = CE_lc(n_tr(1),n_tr(2),p, 1/gamma(kk),D);             
            T2=em_hat_p1*(1+em_hat)*Dg11 + (1+em_hat)^2*Dg12;
             DPc2 = gamma(kk)^2*(1+em_hat)^4*Dg22-2*gamma(kk)*T2+(1+em_hat)^2*Dg11;
             if (DPc2<0)
                epsl_Pc(kk) =  1000;
                cnt=cnt+1;
                continue
             end                 
            G0Pc = (-G0g + 1/n_tr(1)*theta1_hat);
            G1Pc =  (G1g + 1/n_tr(2)*theta1_hat);                      
            epslon_0 = normcdf((G0Pc+cc)/sqrt(DPc2),0,1);  %e0
            epslon_1 = normcdf((G1Pc-cc)/sqrt(DPc2),0,1);  %e1
            epsl_Pc(kk)  =  c0*epslon_0 + c1*epslon_1;   %total error (Proposed consistent)
           end
           [e_min2,indx_opt]=min(epsl_Pc);
           gamma_o=gamma(indx_opt);
            HNL =  V*(D^2/n + gamma_o*eye(nt))^(-2)*D^2*V'/n;