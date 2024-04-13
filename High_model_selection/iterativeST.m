
function [betahat nit] = iterativeST(Y,X,thr,maxit,betahat0,stepsize,nothr);

% iterativeST -- Iterative soft thresholding
%  Usage
%    betahat = iterativeST(Y,X,thr,maxit,betahat0,stepsize);
%  Inputs
%    Y         real vector; observations in sparsity model Y = X*beta+noise 
%    X         real rectangular matrix; model of covariates (design matrix)
%    thr       real number; fixed threshold used in iterations
%    maxit     integer; maximum number of iterations (default is 1000)
%    betahat0  real vector; initial value for betahat (default is zero)
%    stepsize  stepsize in steepest descend direction (default is 1)
%    nothr     integer vector; entries into betahat that are not thresholded
%              (default empty vector)
%  Outputs
%    betahat   estimate of sparse vector beta
%  Description
%    Iterative soft thresholding:
%    betahat_new = ST(betahat_old + stepsize*X'*(Y-X*betahat_old),thr)
%  Notes
%    (1)
%    The default routine assumes matrix X to have induced 2-norm smaller than 1
%    If norm(X) is larger than 1, then use X/c as input instead of X and set
%    betahat = betahat/c, with c = norm(X) (or larger)
%    (Note that this normalisation is different than that in LARS.m)
%    (2)
%    Local optimisation of stepsize is obtained by setting stepsize=NaN
%  See also:
%    help iterativeHT
%    help iterativeGCVST
%    help iterativeMSEST
%    help LARS

Y = column(Y);

if nargin < 2,
   X = eye(length(Y));
end
[n m] = size(X);
if nargin < 3, thr = sqrt(2*log(m)); end
if nargin < 4, maxit = 1000; end
if nargin < 5, betahat0 = zeros(m,1); end
if nargin < 6, stepsize = 1; end
if nargin < 7, nothr = []; end
if max(size(stepsize))>1, nothr = stepsize; stepsize = 1; end
if max(size(betahat0))==1, stepsize = betahat0; betahat0 = zeros(m,1); end
if max(size(maxit))>1, betahat0 = maxit; maxit = 1000; end

beta0 = ones(m,1);
beta1 = betahat0;

if norm(beta1-beta0) < eps, % i.e., if betahat0 happens to be ones(m,1)
   beta0 = zeros(m,1);
end

lambda = thr;

i = 1;
while norm(beta1-beta0) > eps,
   beta0 = beta1;
   dbeta = X'*(Y-X*beta0); % gradient of norm(Y-X*beta) at beta0
   if isnan(stepsize), % then minimize norm(Y-X*beta0-X*alfa*dbeta)
      alfa = norm(dbeta)^2/norm(X*dbeta)^2;
   else
      alfa = stepsize;
   end
   y1 = beta0 + alfa*dbeta;
   beta1 = ST(y1,lambda);
   if length(nothr)>0, beta1(nothr) = y1(nothr); end
   i = i+1; 
   if i>maxit, 
      beta0 = beta1; 
      warning('iterativeST: maximum number of iterations reached')
   end
end
betahat = beta1;
nit = i;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 

