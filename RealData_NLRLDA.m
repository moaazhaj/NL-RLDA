% This is to simulate some examples of binary LDA classification with
% synthetic data to illustrate the performance of NL-RLDA proposed classifier
% by Maaz Mahadi et. al., "Regularized Linear Discriminant Analysis Using a
% Nonlinear Covariance Matrix Estimator"
% 
% The file is a modified version from  "RegularizedSCM" by 
% Esa Ollila and Elias Raninen. 
%  available at: http://users.spa.aalto.fi/esollila/regscm/ 


clear; clc;
addpath (genpath('functions'))
addpath (genpath('datasets'))

load phonemes12.mat

Nvect = [40:50:800];
K       = max(classlabels);  % #of groups(classes)
p       = size(data,2);      % dimension
NRSIM   = 100;               % NR os simulation runs
errLDA  = zeros(1,NRSIM,length(Nvect));   % LDA 

Yf      = (classlabels == 1:K);
nf      = sum(Yf); 
prior = nf./sum(nf);
train_percent_vect = Nvect/sum(nf);


for  jj=1:length(Nvect)
    train_percentage = Nvect(jj)/sum(nf);
     n               =  round(nf*train_percentage); % vector of training samples
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
ylim([15 25])