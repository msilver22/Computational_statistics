
%studnumber = 12345;
studnumber = 202324;
if exist('typeofX') ~= 1, typeofX = 'random'; end
if exist('realrandom') ~= 1, realrandom = false; end
if exist('dimX') ~= 1, dimX = 0; end
   if ~dimX, m = 1000; n = 200; end % n = number of observations
clear dimX
if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
if realrandom==false, rand('state',studnumber); randn('state',studnumber); end

% Set up simulation
p = degreeofsparsity;

beta = zeros(m,1);
I = find(rand(m,1)<p);
beta(I) = randlaplace(size(I),1/5);
%X = randn(n,m)/40; 
pX = 0.1;
X = zeros(n,m); I = find(rand(n,m)<pX); X(I) = (rand(size(I))-0.5)/10;
%X = X*10;
%X= X/10;
%X = X*3;
X= X/3;
mu = X*beta;
stdev = 1;
Z = randn(n,1);
Y = mu+stdev*Z;

maxit = 1000;
univlambda = sqrt(2*log(n))*stdev;

% iterative ST for several thresholds
beta0 = zeros(size(beta));
nlambda = 50;
lambda = (1:nlambda)/nlambda*univlambda/3;
                              % division by factor 3 after having observed that
                              % min PE error is way below univlambda.
samplePEITST = zeros(1,nlambda);
SSeITST = zeros(1,nlambda);
kappa = zeros(1,nlambda);
for t = 1:nlambda
   betahatITST = iterativeST(Y,X,lambda(t),maxit,beta0);
   samplePEITST(t) = norm(X*betahatITST-mu)^2/n;
   kappa(t) = sum(abs(betahatITST)>eps);
   SSeITST(t) = norm(X*betahatITST-Y)^2/n;
end
CpITST = SSeITST+2*kappa*stdev^2/n-stdev^2;

figure(1)
plot(kappa,samplePEITST,'linewidth',3,'color','r')
hold on
plot(kappa,CpITST,'linewidth',2,'color','b')
hold off
title('PE(\kappa) and Cp(\kappa)')
xlabel('model size k')
legend('sample PE','Cp')

figure(2)
plot(lambda,samplePEITST,'linewidth',3,'color','r')
hold on
plot(lambda,CpITST,'linewidth',2,'color','b')
hold off
title('PE(\lambda) and Cp(\lambda)')
xlabel('\lambda')
legend('sample PE','Cp')


[minsamplePE, o] = min(samplePEITST);
lambdaminPE = lambda(o);
kappaminPE = kappa(o);
[minCp, o] = min(CpITST);
lambdaminCp = lambda(o);
kappaminCp = kappa(o);


maxit = 100;
% all steps in iterative soft-thresholding with minimum sample prediction error
[samplePEminPE, lossminPE, SSEminPE] = PElossSSEitST(maxit,Y,X,lambdaminPE,beta,beta0);
% all steps in iterative soft-thresholding with minimum Cp
[samplePEminCp, lossminCp, SSEminCp] = PElossSSEitST(maxit,Y,X,lambdaminCp,beta,beta0);

figure(3)
plot(samplePEminPE(1:maxit),'.-','color','r')
hold on
plot(samplePEminCp(1:maxit),'.-','color','b')
hold off
title('PE(it) and Cp(it)')
xlabel('iteration steps it')
legend('\lambda min PE','\lambda min Cp')



kappatrue = sum(abs(beta)>eps);
kappamax = min(n,kappa(1));
% when kappa>n, then the selected model leads to rank-deficient (singular)
% matrices in the LARS procedure. When kappa is close to n, rank-deficiency may
% occur already. Best is to keep kappamax well below n.

[betahatLARS, XbetahatLARS, selectionLARS, CpLARS, lambdaLARS SLARS] = ...
                                                  LARS(Y,X,stdev,kappamax+1);

figure(4)
plot(kappa,samplePEITST,'linewidth',3,'color','r')
hold on
plot(kappa,CpITST,'linewidth',2,'color','b')
axy = axis;
plot(1:kappamax+1,CpLARS(3:kappamax+3),'linewidth',2,'color','m')
hold off
axis(axy)
legend('sample PE','Cp','CpLARS')
xlabel('model size k')
title('LARS compared to ITST with X = X/3')


[minsamplePE, o] = min(samplePEITST);
lambdaminPE = lambda(o);
kappaminPE = kappa(o);

[minCp, o] = min(CpITST);
lambdaminCp = lambda(o);
kappaminCp = kappa(o);

[minCpLARS, o] = min(CpLARS(3:kappamax+3));
kappaminCpLARS = o;


%Comparing design matrix plotting mu and Y
figure(5)
plot(mu, 'linewidth',1.5,'Color','g'); 
hold on; 
plot(Y,'linewidth',1,'Color' , 'k'); 
hold off; 
legend('\mu', 'Y');
xlabel('components')


