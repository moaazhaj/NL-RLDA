function [em_hat,em_hat_p1b, theta_G_hat] = CE_D2(n0,n1,p,gama,Sigma)
%clear all

n=n0+n1;
M = p;
N = n;
R = Sigma ;


%gama = 0.01;

z = -1/gama;


  I_N = eye(N);
  I_n = eye(n);
  I_M = eye(M);
  I_p = eye(p);
  
  n=n-2;
  N=N-2;
 % C = (Y*Y'/N);
C = R;



%f_gamaC = 1/n*trace ( C*(gama*C + eye(p))^(-1) ); 

fzp1   = 1/n*trace ( C*(C - z*eye(p))^(-1) ); 
fzp2   = 1/n*trace ( C*(C - z*eye(p))^(-2) );
%fzp3   = 2/n*trace ( C*(C - z*eye(p))^(-3) ); 
%fzp4   = 6/n*trace ( C*(C - z*eye(p))^(-4) ); 


f_gamaC2 = 1/n*trace ( C*(gama*C + eye(p))^(-2) ); 

f_gamaC_prime = - 1/n*trace ( C^2*(gama*C + eye(p))^(-2) );
f_zC_prime = 1/n*trace ( C*(C - z*eye(p))^(-2) ); 

% v =  trace ( Theta*C*(C - z*eye(p))^(-2) ); 
% dv =  2*trace ( Theta*C*(C - z*eye(p))^(-3) ); 
% ddv =  6*trace ( Theta*C*(C - z*eye(p))^(-4) );

em_hat     = real(fzp1/(1-fzp1));
em_hat_p1a = real(f_zC_prime/ (1-fzp1)^2);
em_hat_p1b = real((1+em_hat)^2*fzp2 );
%em_hat_p2  = real( (1+em_hat)^2*fzp3 + 2*em_hat_p1b*(1+em_hat)*fzp2 );
%em_hat_p3  = real( (1+em_hat)^2*fzp4 + 2*(1+em_hat)*em_hat_p1b*fzp3+2*fzp3*(em_hat_p1b + em_hat_p1b*em_hat)+2*fzp2*(em_hat_p2+em_hat_p1b^2+em_hat*em_hat_p2));

xm_hat     =  1/(1+em_hat);
xm_hat_p1  =  -em_hat_p1b/(1+em_hat)^2;
%xm_hat_p2  = ( 2*em_hat_p1b^2 - (1+em_hat)*em_hat_p2 )/(1+em_hat)^3;
%xm_hat_p3  = ( (1+em_hat)^2*(4*em_hat_p2*em_hat_p1b-em_hat_p3-em_hat*em_hat_p3-em_hat_p1b*em_hat_p2) - 3*em_hat_p1b*(1+em_hat)*(2*em_hat_p1b^2- em_hat_p2*(1+em_hat) ) ) / (1+em_hat)^5;
%xm_hat_p3b = ( (1+em_hat)^2*(4*ddem-dddem-em*dddem-dem*ddem)- 3*dem*(1+em)*( 2*dem^2-ddem*(1+em)) ) / (1+em)^5

% epsi_tilde_hat_0 = real ( trace (  Theta*C*(C -z*I_p)^(-1) )/(xm_hat));
% epsi_tilde_hat = real ( trace (  Theta*C*(C -z*I_p)^(-2) )/(xm_hat - z*xm_hat_p1));
% epsi_tilde_hat = real ( trace(Theta*C*(C-z*I_M)^(-2))/(xm_hat - z*xm_hat_p1));
% epsi_tilde_hat_p1 = real ( (2*(xm_hat-z*xm_hat_p1)*trace(Theta*C*(C-z*I_M)^(-3))+ trace(Theta*C*(C-z*I_M)^(-2))*z*xm_hat_p2 )/(xm_hat - z*xm_hat_p1)^2);
% nepsi_tilde_hat_p2 = real( (xm_hat-z*xm_hat_p1)*ddv-z*dv*xm_hat_p2 + v*(z*xm_hat_p3 + xm_hat_p2) + z*xm_hat_p2 * dv);
% depsi_tilde_hat_p2 = real ( -2*z*xm_hat_p2*(xm_hat-z*xm_hat_p1));
% epsi_tilde_hat_p2 = real ( ((xm_hat-z*xm_hat_p1)^2*nepsi_tilde_hat_p2 - ((xm_hat-z*xm_hat_p1)*dv + z*v*xm_hat_p2)*depsi_tilde_hat_p2)/(xm_hat-z*xm_hat_p1)^4);

% epsi_tilde_hat_1 = 1/n* trace ( Theta*C*(C -z*I_p)^(-1)) ;
% epsi_tilde_hat_2 = 1/n* trace ( Theta*C*(C -z*I_p)^(-2)) ;
% epsi_tilde_hat_3 = 1/n* trace ( Theta*C*(C -z*I_p)^(-3)) ;
% epsi_tilde_hat_4 = 1/n* trace ( Theta*C*(C -z*I_p)^(-4)) ;


%a1_hat    =  xm_hat-z*xm_hat_p1;
%a2_hat    =  -xm_hat_p1 + (xm_hat_p1^2/xm_hat -0.5*xm_hat_p2)*z;
%a3_hat    =  xm_hat - 2*xm_hat_p1*z +xm_hat_p1^2/xm_hat*z^2;
%a4_hat    =  -0.5*xm_hat_p2 + xm_hat_p1^2/xm_hat + ( -xm_hat_p1^3/xm_hat^2 + xm_hat_p1*xm_hat_p2/xm_hat - xm_hat_p3/6)*z;
%a5_hat    =  -2*xm_hat_p1 + (-xm_hat_p2+ 4*xm_hat_p1^2/xm_hat)*z + (-2*xm_hat_p1^3/xm_hat^2+xm_hat_p1*xm_hat_p2/xm_hat)*z^2;
%a6_hat    =  xm_hat-3*xm_hat_p1*z+3*xm_hat_p1^2/xm_hat*z^2-xm_hat_p1^3/xm_hat^2*z^3;

epsi_hat_2_In = fzp2/(xm_hat-z*xm_hat_p1);
phi_hat       = real(em_hat_p1b-epsi_hat_2_In)/(em_hat_p1b*(1+em_hat)^(-2));
DE2_hat       = (xm_hat-z*xm_hat_p1)*phi_hat;
theta_G_hat   = n*(em_hat_p1b*(xm_hat-z*xm_hat_p1) -fzp2 )/(em_hat_p1b*(1+em_hat)^(-2));
% fdf= 1/n* trace(Theta_hat*(Y*Y'/N)*((Y*Y'/N)-z*I_M)^(-2))-(1/n0+1/n1)*theta_G_hat
%epsi_tilde_1_hat_Sn = 1/n* (1/xm_hat*trace(C) + z/xm_hat*n*em_hat);
%depsi_tilde_hat_2_In  = real ((xm_hat-z*xm_hat_p1)*fzp3 + fzp2*(z*xm_hat_p2))/(xm_hat-z*xm_hat_p1)^2;
%depsi_tilde_hat_2_In  = fzp3/(xm_hat-z*xm_hat_p1) + fzp2*(z*xm_hat_p2)/(xm_hat-z*xm_hat_p1)^2;
%ddepsi_tilde_hat_2_In = fzp4/(xm_hat-z*xm_hat_p1)+fzp3*z*xm_hat_p2/(xm_hat-z*xm_hat_p1)^2 + (fzp2*(z*xm_hat_p3+xm_hat_p2)+ z*xm_hat_p2*fzp3)/(xm_hat-z*xm_hat_p1)^2 + 2*z^2*xm_hat_p2^2*fzp2/(xm_hat-z*xm_hat_p1)^3;
%epsi_tilde_2_hat_Sn   = phi_hat;
%depsi_tilde_2_hat_Sn    =  2*em_hat_p1b/(1+em_hat)^(-1) - (em_hat_p1b*(1+em_hat)^(-2)*depsi_tilde_hat_2_In-epsi_hat_2_In*(-2*em_hat_p1b^2*(1+em_hat)^(-3)+(1+em_hat)^(-2)*em_hat_p2))/(em_hat_p1b^2*(1+em_hat)^(-4))
%depsi_tilde_2_hat_Sn    =  2*em_hat_p1b/(1+em_hat)^(-1) - depsi_tilde_hat_2_In/(em_hat_p1b*(1+em_hat)^(-2))- 2*epsi_hat_2_In/((1+em_hat)^(-1))+ em_hat_p2*epsi_hat_2_In/(em_hat_p1b^2*(1+em_hat)^(-2));
%(2*(1-z*dxm/xm)*epsi_tilde_3-2*dxm/xm*epsi_tilde_2)/n

%ddepsi_tilde_2_hat_Sn = 2*em_hat_p2/(1+em_hat)^(-1) + 2*em_hat_p1b^2 - ddepsi_tilde_hat_2_In/(em_hat_p1b*(1+em_hat)^(-2)) -2*depsi_tilde_hat_2_In/(1+em_hat)^(-1) + em_hat_p2*depsi_tilde_hat_2_In/(em_hat_p1b^2*(1+em_hat)^(-2))...
 %   -2*depsi_tilde_hat_2_In/(1+em_hat)^(-1) -2*em_hat_p1b*epsi_hat_2_In+ (em_hat_p2*depsi_tilde_hat_2_In+em_hat_p3*epsi_hat_2_In)/(em_hat_p1b^2*(1+em_hat)^(-2)) + 2*em_hat_p2*epsi_hat_2_In/(em_hat_p1b*(1+em_hat)^(-1)) - 2*em_hat_p2^2*epsi_hat_2_In/(em_hat_p1b^3*(1+em_hat)^(-2));

%ddepsi_tilde_2_Sn =(3*(1-z*dxm/xm)*epsi_tilde_4-3*dxm/xm*epsi_tilde_3)/n

%depsi_tilde_3_hat_Sn = (ddepsi_tilde_2_hat_Sn + 2*xm_hat_p1/xm_hat*depsi_tilde_2_hat_Sn+2*epsi_tilde_2_hat_Sn*(xm_hat*xm_hat_p2-xm_hat_p1^2)/xm_hat^2)/(2*(1-z*xm_hat_p1/xm_hat))  + ...
%    (depsi_tilde_2_hat_Sn + 2*xm_hat_p1/xm_hat*epsi_tilde_2_hat_Sn)*2*(xm_hat*(xm_hat_p1+xm_hat_p2*z)-xm_hat_p1^2*z)/(4*xm_hat^2*(1-z*xm_hat_p1/xm_hat)^2);

 

%epsi_tilde_3_hat_Sn = (depsi_tilde_2_hat_Sn + 2*xm_hat_p1/xm_hat*epsi_tilde_2_hat_Sn)/(2*(1-z*xm_hat_p1/xm_hat));

 %epsi_tilde_4_hat_Sn = (depsi_tilde_3_hat_Sn + 3*xm_hat_p1/xm_hat*epsi_tilde_3_hat_Sn)/(3*(1-z*xm_hat_p1/xm_hat));

%theta_D2_hat = ((1/xm_hat*(a2_hat+a4_hat*z))*epsi_tilde_hat_1 + (-a1_hat+a3_hat/xm_hat + (a2_hat/xm_hat-2*a2_hat+a5_hat/xm_hat)*z+(a4_hat/xm_hat-a4_hat)*z^2)*epsi_tilde_hat_2 ...
 %   + ((a3_hat/xm_hat-2*a3_hat+a6_hat/xm_hat)*z+(a5_hat/xm_hat-a5_hat)*z^2)*epsi_tilde_hat_3 + ((a6_hat/xm_hat-a6_hat)*z^2)*epsi_tilde_hat_4);
%theta_D2_hat_Sigma= (1/n0+1/n1)*((1/xm_hat*(a2_hat+a4_hat*z))*epsi_tilde_1_hat_Sn + (-a1_hat+a3_hat/xm_hat + (a2_hat/xm_hat-2*a2_hat+a5_hat/xm_hat)*z+(a4_hat/xm_hat-a4_hat)*z^2)*epsi_tilde_2_hat_Sn ...
 %   + ((a3_hat/xm_hat-2*a3_hat+a6_hat/xm_hat)*z+(a5_hat/xm_hat-a5_hat)*z^2)*epsi_tilde_3_hat_Sn + ((a6_hat/xm_hat-a6_hat)*z^2)*epsi_tilde_4_hat_Sn);






