% This is to simulate some examples of binary LDA classification with
% synthetic data to illustrate the performance of NL-RLDA proposed classifier
% by Maaz Mahadi et. al., "Regularized Linear Discriminant Analysis Using a
% Nonlinear Covariance Matrix Estimator"
% 
% The file is a modified version from  "RegularizedSCM" by 
% Esa Ollila and Elias Raninen 
% available at: http://users.spa.aalto.fi/esollila/regscm/ 
% 
clear; clc;
addpath (genpath('functions'))

K           = 2;            % #of groups(classes)
p           = 50;           % dimensionality 
Ntotal      = 10000;        % Population size 
Nvect       = 25:20:105;    % a vector of the number of training samples
c0          = 0.5;
c1          = 1 - c0;
prior       = [c0 , c1];   % prior probabilities (class0,class1) 
NRSIM       = 100;           % Number of simulation runs
errLDA      = zeros(8,NRSIM,length(Nvect));   % LDA 

% Model 1
 Sigma = toeplitz([1 (0.1)*ones(1,p-1)]); 
 
% Model 2
%     x = 1:p;
%     [X,Y] = meshgrid(x,x);
%     Sigma = 0.9.^(abs(X-Y));

% Model 3
%  Sigma = eye(p);
%  subdiag = [repmat(0.9, 1, 4), repmat(0.3, 1, 5)];
%  for k = 1:9
%       Sigma = Sigma + diag(subdiag(k)*ones(1,p-k),k) + diag(subdiag(k)*ones(1,p-k),-k);
%  end
  
  
 [U, SD, V] = svd(Sigma);
 SD_sqrt = diag(sqrt(diag(SD)));
 SD_inv = diag((diag(1./SD)));
 Sigma_sqrt = U * SD_sqrt * V';
 Sigma_inv = U * SD_inv * V';
  
  nu_sq             = 1;
  mu0                = sqrt(nu_sq/(4*sum(sum(Sigma_inv))))*ones(p,1);  
  mu1                = -1*mu0;
  mu                 = [mu0';mu1'];
  nf                 = round(prior*Ntotal);
  train_percent_vect = Nvect/sum(nf);
  Nindf              = [0 cumsum(nf)];

    for ii=1:K
        Zf = randn(p,nf(ii));
        Xf = [mu(ii,:)'*ones(1,nf(ii))] + Sigma_sqrt*Zf;
        data((Nindf(ii)+1):Nindf(ii+1),:) = Xf';       % data
        classlabels((Nindf(ii)+1):Nindf(ii+1))=ii';    %classlabel 
    end



for  jj=1:length(Nvect)
    train_percentage = Nvect(jj)/sum(nf);
     n               =  round(nf*train_percentage); % vector of # of training samples
    for iter=1:NRSIM
         [jj iter]
        % Training set X and test set Z selection by random sampling
         [yt,Xt,y,Xc,hat_mu] = create_data(data,classlabels,train_percentage,true);
        
        % NL_RLDA                   
         H_NL = NL_RLDA(Xc,hat_mu.',n);
         errLDA(1,iter,jj) = LDA_test_error(Xt,yt, hat_mu,H_NL,prior);         
         
    end
end

avgerror = mean(errLDA,2);
avgerror = avgerror(:,:);

plot(Nvect,100*avgerror(1,:),'Linestyle','-','color','r','Marker','.','LineWidth',2, 'MarkerSize',30) 
grid on 
xlabel('No. of samples'), ylabel('Avg. error (%)'), 
AX = legend({'NL-RLDA'},'Location','northwest','NumColumns',1);
legend('boxoff')
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',20); 
set(gca,'FontSize',24);
ylim([25 40])


