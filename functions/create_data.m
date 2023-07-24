function [yt,Xt,y,X,mu,prior] = create_data(Xo,yo,pt,center_train_data)

if nargin < 4
    center_train_data = false;
end
I       = eye(max(yo));   
G       = max(yo);
Yf      = logical(I(yo, :));

% Yf is an n x G  indicator matrix => ind(i,j) = 1, i belong to group j
[Nf,p]  = size(Xo);     % full sample length and dimension
nf      = sum(Yf);      % nf is full sample length in classes
n = round(nf*pt);   % sample lengths in training sets

if any(n == 0) 
    error('ERROR, not enough training data!!!'); 
end

N    = sum(n);  
Nind = [0 cumsum(n)];      
prior = n / N; 	% empirical prior (probability) of classes

%% Generate class label vectors y and yt

ntst    = nf -n;                % sample lengths in test sets
Nt      = sum(ntst);
Nind2   = [0 cumsum(ntst)];      
 
yt = zeros(Nt,1);
y = zeros(N,1);
for ii=1:G
     yt((Nind2(ii)+1):Nind2(ii+1))=ii; % class labels for test set
     y((Nind(ii)+1):Nind(ii+1))=ii; % class labels of training set
end      

%% Generate training data set X and test set Xt

X  = zeros(N,p); 
Xt = zeros(Nt,p); 
mu = zeros(p,G); 
onetoNf = 1:Nf;
 
for ii=1:G
     i_tr  = randsample(onetoNf(Yf(:,ii)),n(ii)); % indices (training)
     i_tst = setdiff(onetoNf(Yf(:,ii)),i_tr);     % indices (testing)
     
     %% Create the training set
     % Note: replace nan's in the training set with the mean values 
     Xii = Xo(i_tr,:);
     ind = isnan(Xii);
     mu(:,ii)  = nanmean(Xii);              % class means
     tmp = repmat(mu(:,ii).',numel(i_tr),1);  
     Xii(ind) = tmp(ind);
     if ~center_train_data 
        X((Nind(ii)+1):Nind(ii+1),:) = Xii;
     else
        X((Nind(ii)+1):Nind(ii+1),:) = Xii - repmat(mu(:,ii).',n(ii),1); 
     end
     
     %% Create test data 
     % replace nan's with corresponding mean values in the test set
     Xii = Xo(i_tst,:);
     ind = isnan(Xii);
     muii  = nanmean(Xii);              % class means
     tmp = repmat(muii,numel(i_tst),1);  
     Xii(ind) = tmp(ind);
     Xt((Nind2(ii)+1):Nind2(ii+1),:) = Xii;% test data set
     %% 
end      
 
 
    