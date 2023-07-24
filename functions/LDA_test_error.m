function  [err_lda] = LDA_test_error(Xt,yt,M,invC,pr)
% LDA_TEST_ERROR computes the test error of the linear discriminant rule
% given the test set and labels, class means, and inverse of the pooled
% sample covariance matrix (SCM)
%
% inputs: 
%     Xt       test data set of size Nt x p. 
%     yt       class labels of test data set (a vector of sized Nt x 1)
%     M        matrix with class means of training data set (of size p x K)
%               as column vectors 
%     invC     p x p matrix of the inverse of the regularized poooled SCM
%
% optional inputs:
%     pr       optional class a priori probabilities (a vector of size 1xK)
%              Elements of pr are in [0,1] and should sum to 1. Defaul value 
%              is pr = (1/K, ..., 1/K). 
%
% Outputs: 
%     err_lda  the misclassification rate
%
% toolbox: RegularizedSCM ({esa.ollila,elias.raninen}@aalto.fi)
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = size(M,2);
Nt = size(Xt,1);

if  nargin < 5 
    pr = 1/K * ones(1,K); % equal prior by default
end

B = invC*M;
const = 0.5*sum(M.*B) - log(pr); % 1xK vector of constant in LDA disc. fnc

%% Predict the labels of training and test set of data
dlda = -Xt*B + repmat(const,Nt,1);% scores
[~,indx_lda] = min(dlda,[],2);
err_lda = sum(indx_lda ~= yt)/Nt;
