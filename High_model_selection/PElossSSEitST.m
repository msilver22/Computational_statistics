
function [PE loss SSE] = PElossSSEitST(nit,Y,X,lambda,beta,beta0)

% PElossSSEitST -- sample PE, loss and SSE as a function of the iteration steps
%  Usage
%    [PE loss Cp] = makeplotPElossCpITST(nit,Y,X,lambda,beta,beta0);
%  Inputs
%    nit       integer; maximum number of iterations
%    Y         real vector; observations in sparsity model Y = X*beta+noise 
%    X         real rectangular matrix; model of covariates (design matrix)
%    lambda    real number; fixed threshold used in iterations
%    beta      real vector; true value
%    beta0     real vector; initial value for betahat
%  Outputs
%    PE        vector of length nit; sample prediction errors
%    loss      vector of length nit
%    SSE       vector of length nit 


SSE = zeros(1,nit);
PE = zeros(1,nit);
loss = zeros(1,nit);
mu = X*beta;
[n m] = size(X);
warning off
for it = 1:nit,
   beta1 = iterativeST(Y,X,lambda,1,beta0);
   loss(it) = norm(beta1-beta)^2/m;
   PE(it) = norm(X*beta1-mu)^2/n;
   SSE(it) = norm(X*beta1-Y)^2/n;
   beta0 = beta1;
end
warning on
