function [HL,gamma_o] = nlrlda_lc(X,mu,n_tr)
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
    N1 = n_tr(1)+n_tr(2);

    cc=log(c1/c0);    
     p = length(S_hat);
   % gg          = -10:0.25:-5;
    gg=linspace(-10,10,21);
    gamma = 1e5.^(gg./10);
    epsl = zeros(1,length(gamma));
            for kk = 1:length(gamma)
                H1 = V*( (gamma(kk)*D^2/(nt-2) + eye(nt))^(-1) - eye(nt))*V' + eye(p);
                %trace_H = trace((gamma(kk)*D^2/(nt-2) + eye(nt))^(-1))-n+p;
                G0 = (mu_hat(:,1) - 1/2*m_plus)'*H1*m_minus;
                G1 = (mu_hat(:,2) - 1/2*m_plus)'*H1*m_minus;
                D0 = m_minus'*H1*S_hat*H1*m_minus;
                del_hat = (p - trace(H1))/(gamma(kk)*(N1-2-p + trace(H1)));
                T1= (1+gamma(kk)*del_hat)^2*D0;
             if (T1<0)
                epsl(kk) =  1000;
                continue
             end 
                eps_0 = normcdf((-G0+(N1-2)*del_hat/n_tr(1)+cc)/sqrt(T1),0,1);  %e0
                eps_1 = normcdf((G1+(N1-2)*del_hat/n_tr(2)-cc)/sqrt(T1),0,1);  %e1
              %  epsl(iter,ii)  =  c0*eps_0 + c1*eps_1;   %total error
                epsl(kk) = c0*eps_0 + c1*eps_1;
            end
            [e_min,indx_opt]=min(epsl);
            gamma_o=gamma(indx_opt);
            HL=V*( (gamma_o*D^2/(nt-2) + eye(nt))^(-1) - eye(nt))*V' + eye(p);